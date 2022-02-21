
// Urban Heat Adaptation Problem, MOEA 
// Author: Rui Shi et al. 2021

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include "./Borg-1.9/borgms.h"
#include <time.h>
#include <algorithm>
#include <typeinfo>
#include <mpi.h>

#define _CRT_SECURE_NO_WARNINGS

// ===================== Define Parameters ================================== // 

using namespace std;

int nvars = 110;                   // 110 decision variables, 10 for each 11 districts 
int nobjs = 5;                     // 5 objectives: life saved, Cost, Equity, Reliability, Carbon Reduction
int nconstrs = 0;                  // 0 constrs

const int Nsamples = 1500;         // 1500 scenarios 
const int Nmodels = 32;            // 32 climate models from NA-Cordex and LOCA
const int Nregions = 11;           // 11 planning districts in Baltimore
const int Nyears = 20;             // 20 years
const int Nsummer = 153;           // 153 days during the summer

double District[Nregions][12];
double ScenMat[Nsamples][17];
double temp[Nmodels][Nsummer][Nyears + 1];
double BWI_temp[Nsummer][5];   // 5 years before simulation start 

int adjacent[11][5] = {
{	1,	2,	3,	10, -1	}, {	0,	9,	6,	2,	7	}, {	3,	0,	1,	7,	4	},
{	10,	5,	4,	2,	0	}, {	7,	2,	3, -1, -1	}, {	10,	3, -1, -1, -1},
{	1,	9,	8, -1, -1	}, {	1,	2,	4, -1, -1	}, {	10,	9,	6, -1, -1	},
{	6,	8,	10,	1, -1	}, {	9,	8,	5,	3,	0	}

};                                  // adjacent matrix, row i defines districts are adjacent to district i

// Scenario independent variables  Table 3 
double tree_area = 30;                                                        // 30 sq.meter per tree 
double FC_tree; double MC_tree = 40;                                          // Tree: FC depends on time; MC: $40/tree/yr
double FC_CP = 2.5; double MC_CP = 0.35;								      // CP: $2.5/m2, $0.35/m2/yr
double FC_CR = 22;  double MC_CR = 0.2;								       	  // CR: $22/m2, $0.2/m2/yr
double OperCost_CC = 2; double Cost_CC = 10000;						     	  // CC: $2/person/day; option cost $10000/center
double NMR65 = 0.000223; double NMR64 = 0.0000111; double NMR = 0.0000836;    // mortality on a non-heat wave day 2.23e-4 1.11e-5
double CCcap = 100;														      // cooling center capacity 
double gamma_factor = 0.0218;											      // carbon reduction (ton) from 1 Tree / Year (48 lb /tree/year)

// Scenario dependent variables Table 2
int ClimModel;                                    // Scenario[][0]					   Climate Model ID
int HWdef;                                        // Scenario[][1]                     Heat Wave Definition 
double TR_tree, TR_CP, TR_CR;                     // Scenario[][2], [3], [4]           Temperature Reduction with Tree, CP, CR
double CCuse, p_vage65, p_vpoor;                  // Scenario[][5], [6], [7]		   CC usage, % visitor>65, % visitor poor
double alpha, beta;                               // Scneario[][8], [9]                % impervious to tree, spillover effect
double APG, AAR, APR, DisR;                       // Scneario[][10], [11], [12], [13]  Annual Population Growth, Annual Aging, Annual Poverty, Real Discount
double HWRR_poor, HWRR_age65, HWRR_age64;         // Scneario[][14], [15], [16]        HWRR for poor, age>65, age<65;

double Tree_cap[Nregions];                        // Available lands for Tree in each region
double CP_cap[Nregions];						  // Available streets for Tree in each region
double CR_cap[Nregions];						  // Available roofs for Tree in each region

// ===================== Self-Defined Function ================================== // 

// Log policy function // hw_state # of hw days in t-1, a,b,d decision variables, LU_cap: total available lands, X: current cumulative lands
double LogPolicy(double hw_state, double a, double b, double d, double LU_cap, double X) {
	double Y = 0;
	Y = LU_cap * max( ( (1 + d) / (1 + exp(-a * (hw_state - b) )) - d), 0.00);
	Y = min(Y, (LU_cap - X));
	return Y;
}

// Temperature reduction functions  
// Tree: non-linear effect:  x, current tree canopy; y, tree canopy at t0; TR_tree, unit tempreature reduction of Tree 
double TRTreeFun(double x, double y) {
	double z = TR_tree * (x * x - y * y);
	return z;
}
// CP: x, cumulative added cool pavement; TR_CP, unit tempreature reduction of CP
double TRCPFun(double x) {
	double z = TR_CP * x;
	return z;
}
// CR: x, cumulative added cool roof; TR_CR, unit tempreature reduction of CR
double TRCRFun(double x) {
	double z = TR_CR * x;
	return z;
}

// Calculate mean temperature given a time period 
double mean_temp(double tmax[Nsummer], int left, int right) {
	double sum = 0.000; double temp_mean = 0.000;
	for (int i = left; i < (right + 1); i++) {
		sum = sum + tmax[i];
	}
	temp_mean = sum / (double)(right - left + 1);
	return(temp_mean);
}

// Heat wave counting function 
  // by AB or Huang 
int hw_count1(double tmax[Nsummer], double temp_threshold, int dur_threshold) {
	int dur_hw = 0;
	int j = 0;

	while (j < (Nsummer-dur_threshold+1)) {
		int z = 0;
		while ( (tmax[j] >= temp_threshold) && (j<Nsummer) ) {
			z = z + 1;
			j = j + 1;
		}

		if (z >= dur_threshold) {
			dur_hw = dur_hw + z;
		}

		j = j + 1;
	}

	return (dur_hw);
}

// heat wave definition by Peng et al. 2011
int hw_count2(double tmax[Nsummer], double temp_threshold1, double temp_threshold2, int dur_threshold) {

	int hw_dur = 0;
	int date[Nsummer] = { 0 }; int j = 0; int flag = 0; int hw_num = 0;

	for (int i = 0; i < Nsummer; i++) {
		if (tmax[i] >= temp_threshold1) {
			date[j] = i;
			j = j + 1;
		}
	}
		
	if (j < dur_threshold) {
		return(0);    // less than 3 days over T1
	}

	int firstday = 0;

	for (int i = 0; i < (j - 2); i++) {
		if (((date[i + 1] - date[i]) == 1) && ((date[i + 2] - date[i + 1]) == 1)) {
			firstday = i;
			break;
		}
		else if (i == (j - 3)) {
			return(0);   // no 3 consecutive days over T1
		}

	}

	int left[Nsummer] = { 0 }, right[Nsummer] = { 0 };
	int left_temp, right_temp = 0; int left_wall = 1;
	int i = firstday;
	while (i < (j - 2)) {
				
		if (((date[i + 1] - date[i]) == 1) && ((date[i + 2] - date[i + 1]) == 1)) {

			hw_num = hw_num + 1;
			left[hw_num] = date[i]; right[hw_num] = date[i + 2];

			while ((tmax[left[hw_num]] > temp_threshold2) && (tmax[right[hw_num]] > temp_threshold2) && (mean_temp(tmax, left[hw_num], right[hw_num]) > temp_threshold1)) {

				left_temp = left[hw_num] - 1; right_temp = right[hw_num] + 1;

				if ((left_temp >= 0) && (tmax[left_temp] > temp_threshold2) && (mean_temp(tmax, left_temp, right[hw_num]) > temp_threshold1)) {
					left[hw_num] = left_temp;
				}

				if ((right_temp <= 152) && (tmax[right_temp] > temp_threshold2) && (mean_temp(tmax, left[hw_num], right_temp) > temp_threshold1)) {
					right[hw_num] = right_temp;
				}

				if ((right_temp != right[hw_num]) && (left_temp != left[hw_num])) {
					break;  // break if no updates 
				}
			}

			if (hw_num >= 2) {
				int t = hw_num - 1;
				while ((left[hw_num] < right[t]) && (t >= left_wall)) {
					hw_dur = hw_dur - (right[t] + 1 - left[t]);   // count the heat wave duration
					t = t - 1;
				}

				if (t < (hw_num - 1)) {
					left_wall = hw_num;
				}
			}

			hw_dur = hw_dur + (right[hw_num] + 1 - left[hw_num]);   // count the heat wave duration

			if (right[hw_num] >= date[j - 1]) {
				return (hw_dur);
			}
			else {
				flag = 0;
				for (int k = 0; k < (j - 2); k++) {      // identify new searching location   
					if (date[k] > right[hw_num]) {
						i = k;
						flag = 1;
						break;
					}
				}

				if (flag == 0) {
					return (hw_dur);
				}

			}

		}
		else {
			i = i + 1;
		}
	}
	return (hw_dur);
}

// sample max
double samplemax(double v[Nsamples]) {
	double T = -9999999999;
	for (int i = 0; i < Nsamples; i++) {
		if (v[i] > T) {
			T = v[i];
		}
	}
	return (T);
}

void CityHEAT(double* vars, double* objs, double* constrs) {

	// Get decision variables from Borg
	double DVmat[11][10] = { 0 };
	for (int r = 0; r < Nregions; r++) {
		DVmat[r][0] = ceil(vars[r * 10]);
		for (int j = 1; j < 10; j++) {
			DVmat[r][j] = vars[r * 10 + j];
		}
	}

	// objective functions 
	double TEM[Nsamples] = { 0 }; 	   		// mortality with adaptation
	double TEM0[Nsamples] = { 0 }; 			// mortliaty without adaptation 
	double REM[Nsamples] = { 0 };  			// mortality reduction
	double Cost[Nsamples] = { 0 }; 			// total cost 
	double Equity[Nsamples] = { 0 };		// equity threshold violation 
	double Reliability[Nsamples] = { 0 };   // relability threshold violation
	double C_Reduce[Nsamples] = { 0 };		// CO2 reduction (ton)

	// intermediate variables 
	double DEM[Nregions][Nyears + 1] = { 0 };
	double DEM0[Nregions][Nyears + 1] = { 0 };
	double OperationCost_CC = 0;
	double temp_r_0[Nsummer] = { 0 };
	double temp_r_t[Nsummer] = { 0 };
	double EM_s[Nregions][Nyears + 1] = { 0 };
	double EM0_s[Nregions][Nyears + 1] = { 0 };
	double Cost_s[Nregions][Nyears + 1] = { 0 };
	int HWD[Nregions][Nyears + 5] = { 0 };
	int HWD0[Nregions][Nyears + 5] = { 0 };
	double HW_state = 0; 

	double Y_Tree[Nregions][Nyears + 1] = { 0 }; double CumTree[Nregions][Nyears + 1] = { 0 }; double CanopyTree[Nregions][Nyears + 1] = { 0 };
	double Y_CP[Nregions][Nyears + 1] = { 0 }; double CumCP[Nregions][Nyears + 1] = { 0 };
	double Y_CR[Nregions][Nyears + 1] = { 0 }; double CumCR[Nregions][Nyears + 1] = { 0 };

	double temp_miti[Nregions][Nyears + 1] = { 0 };
	double temp_miti_spill[Nregions][Nyears + 1] = { 0 };
	double temp_miti1[Nregions][Nyears + 1] = { 0 };

	// starting simulation 
	for (int s = 0; s < Nsamples; s++) {

		// Pass the value of parameters from the ScenMat
		ClimModel = (int)ScenMat[s][0]; HWdef = (int)ScenMat[s][1];
		TR_tree = ScenMat[s][2]; TR_CP = ScenMat[s][3]; TR_CR = ScenMat[s][4];
		CCuse = ScenMat[s][5]; p_vage65 = ScenMat[s][6]; p_vpoor = ScenMat[s][7];
		alpha = ScenMat[s][8]; beta = ScenMat[s][9];
		APG = ScenMat[s][10]; AAR = ScenMat[s][11]; APR = ScenMat[s][12]; DisR = ScenMat[s][13];
		HWRR_poor = ScenMat[s][14]; HWRR_age65 = ScenMat[s][15]; HWRR_age64 = ScenMat[s][16];
		unsigned int flag[Nregions] = { 0 };

		// build up HWD matrix for first 5 years
		for (int r = 0; r < Nregions; r++) {
			
			for (int tt = 0; tt < 5; tt++) {

				for (int j = 0; j < Nsummer; j++) {
					temp_r_0[j] = BWI_temp[j][tt] + District[r][6];                     // District[r][6] = dT[r]
				}

				// which heaw wave function to use 
				if (HWdef == 1) {										 // ScenMat[s][2] index for heat wave definition
					HWD[r][tt] = hw_count1(temp_r_0, 35, 3);				 //Huang 2010, Daily maximum temperature > 35 degree C for at least 3 consecutive days
					HWD0[r][tt] = HWD[r][tt];
				}
				else if (HWdef == 2) {
					HWD[r][tt] = hw_count1(temp_r_0, 35.00, 2);		     // AB 2011 Daily maximum temperature > 95th percentile for >= 2 consecutive days
					HWD0[r][tt] = HWD[r][tt];
				}
				else {
					HWD[r][tt] = hw_count2(temp_r_0, 36.10, 32.20, 3);    // Peng  Double threshold T1 = 97.5% and T2 = 81%
					HWD0[r][tt] = HWD[r][tt];
				}
			}

			Cost[s] = Cost[s] + DVmat[r][0] * Cost_CC;   // upfront cost of cooling center 

			CumTree[r][0] = District[r][7]; CanopyTree[r][0] = District[r][7]; 
			Tree_cap[r] = CumTree[r][0] + District[r][10] + alpha * District[r][11];
			CP_cap[r] = District[r][9]; CR_cap[r] = District[r][8];
		}

		// starting from the year 2020 
		for (int t = 1; t < (Nyears+1); t++) {

			// calculate new-installed tree, CP, and CR, update CumTree, CanopyTree,  CumCP, CumCR, calculate temp_miti without spillover effect 
			for (int r = 0; r < Nregions; r++) {
				// moving average of the previous 5 years 
				HW_state = ((double)(HWD[r][t - 1] + HWD[r][t] + HWD[r][t+1] + HWD[r][t+2] + HWD[r][t+3])) / 5;
				Y_Tree[r][t] = LogPolicy(HW_state, DVmat[r][1], DVmat[r][2], DVmat[r][3], Tree_cap[r], CumTree[r][t - 1]);
				Y_CP[r][t] = LogPolicy(HW_state, DVmat[r][4], DVmat[r][5], DVmat[r][6], CP_cap[r], CumCP[r][t - 1]);
				Y_CR[r][t] = LogPolicy(HW_state, DVmat[r][7], DVmat[r][8], DVmat[r][9], CR_cap[r], CumCR[r][t - 1]);
				CumTree[r][t] = CumTree[r][t - 1] + Y_Tree[r][t];
				CanopyTree[r][t] = CanopyTree[r][t - 1] + (CumTree[r][t] - CumTree[r][0]) / 20;
				CumCP[r][t] = CumCP[r][t - 1] + Y_CP[r][t];
				CumCR[r][t] = CumCR[r][t - 1] + Y_CR[r][t];
				temp_miti[r][t] = TRTreeFun(CanopyTree[r][t], CanopyTree[r][0]) + TRCPFun(CumCP[r][t]) + TRCRFun(CumCR[r][t]); // temp mitigation in year t without spillover effect
			}

			// calulate temp_red = temp_miti + spillover effect
			for (int r = 0; r < Nregions; r++) {

				int j = 0;
				while ((adjacent[r][j] >= 0) && (j < 5)) {
					temp_miti_spill[r][t] = temp_miti_spill[r][t] + beta * District[adjacent[r][j]][0] / District[r][0] * temp_miti[adjacent[r][j]][t];  // spillover effect
					j = j + 1;
				}
				temp_miti1[r][t] = temp_miti_spill[r][t] + temp_miti[r][t];  // spillover effect + original reduction 
				temp_miti_spill[r][t] = 0;                                   // clear temp_red for next iteration

				for (int j = 0; j < Nsummer; j++) {
					temp_r_t[j] = temp[ClimModel][j][t] + District[r][6] - temp_miti1[r][t];  // District[r][6] = dT[r]
					temp_r_0[j] = temp[ClimModel][j][t] + District[r][6];                     // District[r][6] = dT[r]
				}

				// which heaw wave function to use 
				if (HWdef == 1) {										 // ScenMat[s][2] index for heat wave definition
					HWD[r][(t+4)]  = hw_count1(temp_r_t, 35, 3);				 //Huang 2010, Daily maximum temperature > 35 degree C for at least 3 consecutive days
					HWD0[r][(t + 4)] = hw_count1(temp_r_0, 35, 3);
				}
				else if (HWdef == 2) {
					HWD[r][(t + 4)]  = hw_count1(temp_r_t, 35.00, 2);		     // AB 2011 Daily maximum temperature > 95th percentile for >= 2 consecutive days
					HWD0[r][(t + 4)] = hw_count1(temp_r_0, 35.00, 2);
				}
				else {
					HWD[r][(t + 4)]  = hw_count2(temp_r_t, 36.10, 32.20, 3);    // Peng  Double threshold T1 = 97.5% and T2 = 81%
					HWD0[r][(t + 4)] = hw_count2(temp_r_0, 36.10, 32.20, 3);
				}

				OperationCost_CC = 0;

				if (HWD0[r][(t + 4)] > 0) {

					// demographics 
					double pop = District[r][1] * pow((1 + APG), t);                       // population at time t,  District[r][1] population at year 0
					double p_age65 = min(1.00, (District[r][3] + AAR * t));				   // % over 65 at time t,   District[r][3] % over 65 at year 0
					double p_poor = min(1.00, (District[r][5] + APR * t));				   // % poor at time t,  District[r][5] % poor at year 0

					double pop_age65poor = pop * p_age65 * p_poor;                         // population over 65 & poor at time t 
					double pop_age65rich = pop * p_age65 * (1 - p_poor);				   // population over 65 & rich at time t
					double pop_age64poor = pop * (1 - p_age65) * p_poor;				   // population under 65 & poor at time t
					double pop_age64rich = pop * (1 - p_age65) * (1 - p_poor);			   // population under 65 & rich at time t

					double HWMR_age65poor = NMR65 * HWRR_age65 * HWRR_poor;				   // heat wave mortality of people over 65 & poor
					double HWMR_age65rich = NMR65 * HWRR_age65;							   // heat wave mortality of people over 65 & poor
					double HWMR_age64poor = NMR64 * HWRR_age64 * HWRR_poor;				   // heat wave mortality of people over 65 & poor
					double HWMR_age64rich = NMR64 * HWRR_age64;						       // heat wave mortality of people over 65 & poor

					DEM0[r][t] = pop_age65poor * HWMR_age65poor + pop_age64poor * HWMR_age64poor +
						pop_age65rich * HWMR_age65rich + pop_age64rich * HWMR_age64rich;
					
					if (HWD[r][(t + 4)] > 0) {
						// calculate cooling center users 
						double	user_age65poor = min(pop_age65poor, DVmat[r][0] * CCcap * CCuse * p_vage65 * p_vpoor);
						double	user_age65rich = min(pop_age65rich, DVmat[r][0] * CCcap * CCuse * p_vage65 * (1 - p_vpoor));
						double 	user_age64poor = min(pop_age64poor, DVmat[r][0] * CCcap * CCuse * (1 - p_vage65) * p_vpoor);
						double	user_age64rich = min(pop_age64rich, DVmat[r][0] * CCcap * CCuse * (1 - p_vage65) * (1 - p_vpoor));
						double  users = user_age65poor + user_age65rich + user_age64poor + user_age64rich;

						// calculate mortality of each group
						DEM[r][t] = (pop_age65poor - user_age65poor) * HWMR_age65poor + (pop_age65rich - user_age65rich) * HWMR_age65rich +
							(pop_age64poor - user_age64poor) * HWMR_age64poor + (pop_age64rich - user_age64rich) * HWMR_age64rich;

						OperationCost_CC = OperCost_CC * users * HWD[r][(t + 4)];
					}
				}

				EM0_s[r][t] = DEM0[r][t] * (double)HWD0[r][(t + 4)];
				EM_s[r][t] = DEM[r][t] * (double)HWD[r][(t + 4)];

				if (flag[r] == 0) {
					if ((CumTree[r][t] - CumTree[r][0]) > District[r][10]) {

						FC_tree = (400 * (District[r][10] + CumTree[r][0] - CumTree[r][t - 1]) +
							600 * (CumTree[r][t] - District[r][10] - CumTree[r][0])) / Y_Tree[r][t];
						flag[r] = 1;
					}
					else {
						FC_tree = 400;
					}
				}
				else {
					FC_tree = 600; // $600 per tree
				}

				double Ytree_num = District[r][0] * Y_Tree[r][t] * 2589988 / tree_area;
				double Cumtree_num = District[r][0] * (CumTree[r][t] - CumTree[r][0]) * 2589988 / tree_area;
				double YCP_sqmeter = District[r][0] * Y_CP[r][t] * 2589988;
				double CumCP_sqmeter = District[r][0] * CumCP[r][t] * 2589988;
				double YCR_sqmeter = District[r][0] * Y_CR[r][t] * 2589988;
				double CumCR_sqmeter = District[r][0] * CumCR[r][t] * 2589988;

				Cost_s[r][t] = (FC_tree * Ytree_num + MC_tree * Cumtree_num + FC_CP * YCP_sqmeter + MC_CP * CumCP_sqmeter +
					FC_CR * YCR_sqmeter + MC_CR * CumCR_sqmeter + OperationCost_CC) / pow((1 + DisR), t);

			}

		}

		// summarize scenario S 
		double max_EM = -99999999;                           // max one-year expected mortality in a given scenario
		double max_EMdiff = -99999999;                       // max one-year mortality rate difference in a given scenario
		for (int t = 1; t < (Nyears + 1); t++) {

			double EM_t_max = -99999999; double EM_t_min = 99999999;
			double EM_t = 0;

			for (int r = 0; r < Nregions; r++) {
				TEM[s] = TEM[s] + EM_s[r][t];
				TEM0[s] = TEM0[s] + EM0_s[r][t];
				EM_t = EM_t + EM_s[r][t];
				C_Reduce[s] = C_Reduce[s] + gamma_factor * (CanopyTree[r][t] - CanopyTree[r][0]) * District[r][0] * 2589988 / tree_area;
				Cost[s] = Cost[s] + Cost_s[r][t];

				double EM_r_pect = EM_s[r][t] / (District[r][1] * pow((1 + APG), t));   // mortality differences 
				if (EM_r_pect > EM_t_max) {
					EM_t_max = EM_r_pect;
				}

				if (EM_r_pect < EM_t_min) {
					EM_t_min = EM_r_pect;
				}
			}

			if ((EM_t_max - EM_t_min) > max_EMdiff) {
				max_EMdiff = (EM_t_max - EM_t_min);
			} // update the max_EMdiff if the difference is larger than the previous one 

			if (EM_t > max_EM) {
				max_EM = EM_t;
			} // update the max_EM if the difference is larger than the previous one 

		}
				
		REM[s] = TEM0[s] - TEM[s];      // annual mortaity reduction 
		Equity[s] = max_EMdiff;         // maximum of one-year mortality difference  
		Reliability[s] = max_EM;        // maximum of one-year mortality 

	}

	// summarize objectives for Borg Optimization 
	for (int s = 0; s < Nsamples; s++) {
		objs[0] = objs[0] + REM[s];
		objs[1] = objs[1] + Cost[s];
		objs[2] = objs[2] + Equity[s];
		objs[3] = objs[3] + Reliability[s];
		objs[4] = objs[4] + C_Reduce[s];
	}

	objs[0] = -objs[0] / Nsamples;   // max
	objs[1] =  objs[1] / Nsamples;
	objs[2] =  objs[2] / Nsamples;
	objs[3] =  objs[3] / Nsamples;
	objs[4] = -objs[4] / Nsamples;  // max 
}

// ===================== Main Function ================================== // 

int main(int argc, char* argv[]) {

	// Main FUN. Part I: INPUT file read

	/*District Mat: 0.Area (sq.mi); 1.Population; 2.Age>=65; 3.Age>=65(%)
	4.Poor; 5.Poor(%); 6.dT; 7.TreePect; 8.RoofPect; 9.RoadPect; 10.TreeAvailPect
	11. ImpSuf NonBuild Pect */

	ifstream infile;
	infile.open("../Data/District.txt");
	if (!infile) {
		perror("Error");
	}
	else {
		for (int i = 0; i < Nregions; i++) {
			for (int j = 0; j < 12; j++) {
				infile >> District[i][j];
			}
		}
		infile.close();
	}

	/*Scenario Mat 1500 by 17 */
	ifstream infile1;
	infile1.open("../Data/Scenarios.txt");
	if (!infile1) {
		perror("Error");
	}
	else {
		for (int i = 0; i < Nsamples; i++) {
			for (int j = 0; j < 17; j++) {
				infile1 >> ScenMat[i][j];
			}
		}
		infile1.close();
	}

	// read the cliamte model temperature series 
	ifstream infile2;
	infile2.open("../Data/temp20202039.txt");
		if (!infile2) {
		perror("Error");
	}
	else {
		for (int i = 0; i < Nmodels; i++)
		{
			for (int k = 0; k < (Nyears + 1); k++) {
				for (int j = 0; j < Nsummer; j++) {
					infile2 >> temp[i][j][k];
				}
			}
		}
		infile2.close();
	}

    // read the heat wave days matrixs for the pervious 5 years.  
	ifstream infile3;
	infile3.open("../Data/BWI_temp.txt");
		if (!infile3) {
		perror("Error");
	}
	else {
		for (int i = 0; i < Nsummer; i++)
		{
			for (int j = 0; j < 5; j++) {
				infile3 >> BWI_temp[i][j];
			}
				
		}
		infile3.close();
	}

	// Main FUN. Part II: BORG setup

	// setting random seed 
	unsigned int seed = atoi(argv[1]);
	srand(seed);
	int NFE = atoi(argv[2]);

	// interface with Borg-MS
	BORG_Algorithm_ms_startup(&argc, &argv);
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE / 200);     // outupt every 200 iterations

	// Define the problem.
	BORG_Problem problem = BORG_Problem_create(nvars, nobjs, 0, CityHEAT);

	// Set the epsilon values used by the Borg MOEA.  Epsilons define the
	// Set all the parameter bounds and epsilons
	for (int j = 0; j < (nvars / 10); j++) {
		BORG_Problem_set_bounds(problem, j * 10,      0.0, 10.0);
		BORG_Problem_set_bounds(problem, j * 10 + 1,  0.0, 10.0);        // tree: a
		BORG_Problem_set_bounds(problem, j * 10 + 2,  0.0, 200.0);       // tree: b
		BORG_Problem_set_bounds(problem, j * 10 + 3, -1.0, 1.0);         // tree: d  
		BORG_Problem_set_bounds(problem, j * 10 + 4,  0.0, 10.0);        //   CP: a
		BORG_Problem_set_bounds(problem, j * 10 + 5,  0.0, 200.0);       //   CP: b
		BORG_Problem_set_bounds(problem, j * 10 + 6, -1.0, 1.0);         //   CP: d
		BORG_Problem_set_bounds(problem, j * 10 + 7,  0.0, 10.0);        //   CR: a
		BORG_Problem_set_bounds(problem, j * 10 + 8,  0.0, 200.0);       //   CR: b
		BORG_Problem_set_bounds(problem, j * 10 + 9, -1.0, 1.0);		 //   CR: d
	}

	BORG_Problem_set_epsilon(problem, 0, 10);       // average toatl mortality
	BORG_Problem_set_epsilon(problem, 1, 1e7);      // average total cost
	BORG_Problem_set_epsilon(problem, 2, 1e-5);     // largest mortality rate differences
	BORG_Problem_set_epsilon(problem, 3, 10);       // worst-year mortality
	BORG_Problem_set_epsilon(problem, 4, 5e4);      // average carbon reduction

	// This is set up to run only one seed at a time 
	char outputFilename[256];
	char runtime[256];
	FILE* outputFile = NULL;
	sprintf(outputFilename, "../Optimization/sets/UrbanHeat_S%d.txt", seed);   // name the set file with UrbanHeat_S%d.set format
	sprintf(runtime, "../Optimization/runtime/UrbanHeat_S%d.runtime", seed);   // name the runtime file with UrbanHeat_S%d.runtime format

	BORG_Algorithm_output_runtime(runtime);

	BORG_Random_seed(seed);
	BORG_Archive result = BORG_Algorithm_ms_run(problem);   // this actually runs the optimization 

	// if this is the master node, print out the final archive 
	if (result != NULL) {
		outputFile = fopen(outputFilename, "w");
		if (!outputFile) {
			BORG_Debug("Unable to open final output file\n");
		}
		BORG_Archive_print(result, outputFile);
		BORG_Archive_destroy(result);
		fclose(outputFile);
	}

	// Shutdown the parallel processes and exit.
	BORG_Algorithm_ms_shutdown();
	BORG_Problem_destroy(problem);
	
	return EXIT_SUCCESS;
}

