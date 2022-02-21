
// Urban Heat Adaptation Problem, MOEA 
// Author: Rui Shi et al.

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <typeinfo>
#include "moeaframework.h"

#define _CRT_SECURE_NO_WARNINGS

// ===================== Define Parameters ================================== // 

using namespace std;

int nvars = 110;                   // 110 decision variables, 10 for each 11 districts 
int nobjs = 11;                    // 11 objectives
int nconstrs = 0;                  // 0 constrs

const int Nsamples = 1500;         // 1500 re-valuation scenarios, change to 3300 for out-of-sample re-evaluation, see line 648
const int Nmodels = 32;            // 32 climate models from NA-Cordex and LOCA
const int Nregions = 11;           // 11 planning districts in Baltimore
const int Nyears = 20;             // 20 years
const int Nsummer = 153;           // 153 days during the summer

double District[Nregions][12];
double ScenMat[Nsamples][17];
double temp[Nmodels][Nsummer][Nyears + 1];
double BWI_temp[Nsummer][5] = { 0 };

int adjacent[11][5] = {
{	1,	2,	3,	10, -1	}, {	0,	9,	6,	2,	7	}, {	3,	0,	1,	7,	4	},
{	10,	5,	4,	2,	0	}, {	7,	2,	3, -1, -1	}, {	10,	3, -1, -1, -1},
{	1,	9,	8, -1, -1	}, {	1,	2,	4, -1, -1	}, {	10,	9,	6, -1, -1	},
{	6,	8,	10,	1, -1	}, {	9,	8,	5,	3,	0	}

};  // adjacent matrix, row i defines districts are adjacent to district i

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

// file to output details
char output_file_name[256];
FILE* outputFile = NULL;

// ===================== Self-Defined Function ================================== // 

// Log policy function 
double LogPolicy(double hw_state, double a, double b, double d, double LU_cap, double X) {
	double Y = 0;
	Y = LU_cap * max( ( (1 + d) / (1 + exp(-a * (hw_state - b) )) - d), 0.00);
	Y = min(Y, (LU_cap - X));
	return Y;
}

// Temperature reduction functions 
double TRTreeFun(double x, double y) {
	double z = TR_tree * (x * x - y * y);
	return z;
}

double TRCPFun(double x) {
	double z = TR_CP * x;
	return z;
}

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

	while (j < Nsummer) {
		int z = 0;
		while (tmax[j] >= temp_threshold) {
			z = z + 1;
			if (j == (Nsummer - 2)) {
				break;
			}
			else {
				j = j + 1;

			}
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

	int firstday =0;

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

// GINI index 
double GINI_index(double EM_pert[Nregions]) {
	double S = 0; double sum_EM_pert = 0; double gini = 0;
	for (int i = 0; i < Nregions; i++) {
		for (int j = 0; j < Nregions; j++) {
			S = S + abs(EM_pert[i] - EM_pert[j]);
		}
	}

	for (int i = 0; i < Nregions; i++) {
		sum_EM_pert = sum_EM_pert + EM_pert[i];
	}

	gini = S / (2 * Nregions * sum_EM_pert);

	return(gini);
}

// GINI index 2
double GINI_index2(double EM_pert[4]) {
	double S = 0; double sum_EM_pert = 0; double gini = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			S = S + abs(EM_pert[i] - EM_pert[j]);
		}
	}

	for (int i = 0; i < 4; i++) {
		sum_EM_pert = sum_EM_pert + EM_pert[i];
	}

	gini = S / (2 * 4 * sum_EM_pert);

	return(gini);
}

// Simulation function: simulate heat wave mortality and different adapataion measures 
void CityHEAT(double* vars, double* objs, double* constrs) {
	
	double DVmat[11][10] = { 0 };

	for (int r = 0; r < Nregions; r++) {
		DVmat[r][0] = ceil(vars[r * 10]);
		for (int j = 1; j < 10; j++) {
			DVmat[r][j] = vars[r * 10 + j];
		}
	}

	// objective functions 
	double TEM=0; 	   	    	// mortality with adaptation
	double TEM0=0; 			    // mortliaty without adaptation 
	double REM=0;  			    // mortality reduction
	double Cost=0; 			    // total cost 
	double Reliability=0;       // relability threshold violation
	double C_Reduce=0;	 	    // CO2 reduction (ton)
	double Equity = 0;   		// equity1 max max_diff 
	double Equity2 = 0;         // equity2 avg max_diff
	double Equity3 = 0;         // equiyt3 max gini 
	double Equity4 = 0;         // equity4 avg gini 

	double OperationCost_CC = 0;	
	double temp_r_0[Nsummer] = { 0 };  
	double temp_r_t[Nsummer] = { 0 };
	double EM_s[Nregions][Nyears + 1] = { 0 }; 
	double EM_pect[Nregions][Nyears + 1] = { 0 };
	double inv_PCap[Nregions][Nyears + 1] = { 0 };
	double EM0_s[Nregions][Nyears + 1] = { 0 };
	double Cost_s[Nregions][Nyears + 1] = { 0 };
	int HWD[Nregions][Nyears + 5] = { 0 }; 
	int HWD0[Nregions][Nyears + 5] = { 0 };
	double DEM[Nregions][Nyears + 1] = { 0 }; 
	double DEM0[Nregions][Nyears + 1] = { 0 };
	double HW_state[Nregions][Nyears + 1] = { 0 };

	double Y_Tree[Nregions][Nyears + 1] = { 0 }; double CumTree[Nregions][Nyears + 1] = { 0 }; double CanopyTree[Nregions][Nyears + 1] = { 0 };
	double Y_CP[Nregions][Nyears + 1] = { 0 }; double CumCP[Nregions][Nyears + 1] = { 0 };
	double Y_CR[Nregions][Nyears + 1] = { 0 }; double CumCR[Nregions][Nyears + 1] = { 0 };
	double Cost_Tree_t[Nregions][Nyears + 1] = { 0 }; double Cost_CP_t[Nregions][Nyears + 1] = { 0 }; double Cost_CR_t[Nregions][Nyears + 1] = { 0 };
	double Cost_CC_t[Nregions][Nyears + 1] = { 0 };

	double temp_miti[Nregions][Nyears + 1] = { 0 };
	double temp_miti_spill[Nregions][Nyears + 1] = { 0 };
	double temp_miti1[Nregions][Nyears + 1] = { 0 };

	double HWMR_age65poor = NMR65 * HWRR_age65 * HWRR_poor;				   // heat wave mortality of people over 65 & poor
	double HWMR_age65rich = NMR65 * HWRR_age65;							   // heat wave mortality of people over 65 & poor
	double HWMR_age64poor = NMR64 * HWRR_age64 * HWRR_poor;				   // heat wave mortality of people over 65 & poor
	double HWMR_age64rich = NMR64 * HWRR_age64;						       // heat wave mortality of people over 65 & poor

	double pop[Nregions][Nyears + 1] = { 0 };                              // population at time t,  District[r][1] population at year 0
	double p_age65[Nregions][Nyears + 1] = { 0 };				           // % over 65 at time t,   District[r][3] % over 65 at year 0
	double p_poor[Nregions][Nyears + 1] = { 0 };				           // % poor at time t,  District[r][5] % poor at year 0

	double pop_age65poor[Nregions][Nyears + 1] = { 0 };                    // population over 65 & poor at time t 
	double pop_age65rich[Nregions][Nyears + 1] = { 0 };				       // population over 65 & rich at time t
	double pop_age64poor[Nregions][Nyears + 1] = { 0 };				       // population under 65 & poor at time t
	double pop_age64rich[Nregions][Nyears + 1] = { 0 };			           // population under 65 & rich at time t

	double EM0_age65poor[Nregions][Nyears + 1] = { 0 };
	double EM0_age64poor[Nregions][Nyears + 1] = { 0 };
	double EM0_age65rich[Nregions][Nyears + 1] = { 0 };
	double EM0_age64rich[Nregions][Nyears + 1] = { 0 };

	double EM_age65poor[Nregions][Nyears + 1] = { 0 };
	double EM_age64poor[Nregions][Nyears + 1] = { 0 };
	double EM_age65rich[Nregions][Nyears + 1] = { 0 };
	double EM_age64rich[Nregions][Nyears + 1] = { 0 };

	// starting simulation 
	unsigned int flag[Nregions] = { 0 };

	// count for the heat wave numbers in year 2019 (t=0)
	for (int r = 0; r < Nregions; r++) {

		// calculate the HWD of the first 5 years
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

		Cost = Cost + DVmat[r][0] * Cost_CC;   // upfront cost of cooling center 
		CumTree[r][0] = District[r][7]; CanopyTree[r][0] = District[r][7]; Tree_cap[r] = CumTree[r][0] + District[r][10] + alpha * District[r][11];
		CP_cap[r] = District[r][9]; CR_cap[r] = District[r][8];

	}

	// caculate the total mortatlity WITH & WITHOUT adaptation

	for (int t = 1; t < (Nyears + 1); t++) {  // starting from 2nd year 2020, using t-1 year's HWD as indication for adaptation  

		// calculate new-installed tree, CP, and CR, update CumTree, CanopyTree,  CumCP, CumCR, calculate temp_miti without spillover effect 
		for (int r = 0; r < Nregions; r++) {
			HW_state[r][t] = ((double)(HWD[r][t - 1] + HWD[r][t] + HWD[r][t + 1] + HWD[r][t + 2] + HWD[r][t + 3])) / 5;
			Y_Tree[r][t] = LogPolicy(HW_state[r][t], DVmat[r][1], DVmat[r][2], DVmat[r][3], Tree_cap[r], CumTree[r][t - 1]);
			Y_CP[r][t] = LogPolicy(HW_state[r][t], DVmat[r][4], DVmat[r][5], DVmat[r][6], CP_cap[r], CumCP[r][t - 1]);
			Y_CR[r][t] = LogPolicy(HW_state[r][t], DVmat[r][7], DVmat[r][8], DVmat[r][9], CR_cap[r], CumCR[r][t - 1]);
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
				HWD[r][(t+4)] = hw_count1(temp_r_t, 35, 3);				 //Huang 2010, Daily maximum temperature > 35 degree C for at least 3 consecutive days
				HWD0[r][(t + 4)] = hw_count1(temp_r_0, 35, 3);
			}
			else if (HWdef == 2) {
				HWD[r][(t + 4)] = hw_count1(temp_r_t, 35.00, 2);		     // AB 2011 Daily maximum temperature > 95th percentile for >= 2 consecutive days
				HWD0[r][(t + 4)] = hw_count1(temp_r_0, 35.00, 2);
			}
			else {
				HWD[r][(t + 4)] = hw_count2(temp_r_t, 36.10, 32.20, 3);    // Peng  Double threshold T1 = 97.5% and T2 = 81%
				HWD0[r][(t + 4)] = hw_count2(temp_r_0, 36.10, 32.20, 3);
			}
			 
			OperationCost_CC = 0;
			// demographics 
			pop[r][t] = District[r][1] * pow((1 + APG), t);                        // population at time t,  District[r][1] population at year 0
			p_age65[r][t] = min(1.00, (District[r][3] + AAR * t));				   // % over 65 at time t,   District[r][3] % over 65 at year 0
			p_poor[r][t] = min(1.00, (District[r][5] + APR * t));				   // % poor at time t,  District[r][5] % poor at year 0

			pop_age65poor[r][t] = pop[r][t] * p_age65[r][t] * p_poor[r][t];                         // population over 65 & poor at time t 
			pop_age65rich[r][t] = pop[r][t] * p_age65[r][t] * (1 - p_poor[r][t]);				   // population over 65 & rich at time t
			pop_age64poor[r][t] = pop[r][t] * (1 - p_age65[r][t]) * p_poor[r][t];				   // population under 65 & poor at time t
			pop_age64rich[r][t] = pop[r][t] * (1 - p_age65[r][t]) * (1 - p_poor[r][t]);			   // population under 65 & rich at time t

			if (HWD0[r][(t + 4)] > 0) {
								
				EM0_age65poor[r][t] = pop_age65poor[r][t] * HWMR_age65poor;
				EM0_age64poor[r][t] = pop_age64poor[r][t] * HWMR_age64poor;
				EM0_age65rich[r][t] = pop_age65rich[r][t] * HWMR_age65rich;
				EM0_age64rich[r][t] = pop_age64rich[r][t] * HWMR_age64rich;

				DEM0[r][t] = EM0_age65poor[r][t] + EM0_age64poor[r][t] + EM0_age65rich[r][t] + EM0_age64rich[r][t];

		        if (HWD[r][(t + 4)] > 0){
		        	// calculate cooling center users 
					double	user_age65poor = min(pop_age65poor[r][t], DVmat[r][0] * CCcap * CCuse * p_vage65 * p_vpoor);
					double	user_age65rich = min(pop_age65rich[r][t], DVmat[r][0] * CCcap * CCuse * p_vage65 * (1 - p_vpoor));
					double 	user_age64poor = min(pop_age64poor[r][t], DVmat[r][0] * CCcap * CCuse * (1 - p_vage65) * p_vpoor);
					double	user_age64rich = min(pop_age64rich[r][t], DVmat[r][0] * CCcap * CCuse * (1 - p_vage65) * (1 - p_vpoor));
					double  users = user_age65poor + user_age65rich + user_age64poor + user_age64rich;

					EM_age65poor[r][t] = (pop_age65poor[r][t] - user_age65poor) * HWMR_age65poor;
					EM_age64poor[r][t] = (pop_age64poor[r][t] - user_age64poor) * HWMR_age64poor;
					EM_age65rich[r][t] = (pop_age65rich[r][t] - user_age65rich) * HWMR_age65rich;
					EM_age64rich[r][t] = (pop_age64rich[r][t] - user_age64rich) * HWMR_age64rich;

					// calculate mortality of each group
					DEM[r][t] = EM_age65poor[r][t] + EM_age64poor[r][t] + EM_age65rich[r][t] + EM_age64rich[r][t];

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
			Cost_Tree_t[r][t] = (FC_tree * Ytree_num + MC_tree * Cumtree_num) / pow((1 + DisR), t);
			Cost_CP_t[r][t] = (FC_CP * YCP_sqmeter + MC_CP * CumCP_sqmeter) / pow((1 + DisR), t);
			Cost_CR_t[r][t] = (FC_CR * YCR_sqmeter + MC_CR * CumCR_sqmeter) / pow((1 + DisR), t);
			Cost_CC_t[r][t] = (OperationCost_CC) / pow((1 + DisR), t);
		}


	}

	// summaeize scenario S 
	double max_EM = -99999999; 
	double max_max_EMdiff = -99999999;     // Equity 	
	double mean_max_EMdiff = 0;            // Equity2
	double max_gini_t = -99999999;         // Equity3 
	double mean_gini_t = 0;                // Equity4 

	for (int t = 1; t < (Nyears + 1); t++) {

		double EM_t_max = -99999999; double EM_t_min = 99999999;
		double EM_t = 0; double EM0_t = 0; double Cost_t = 0; double C_Reduce_t = 0;
		double gini_t = 0; double EM_r_pect[Nregions] = { 0 }; 
		double gini_inv = 0; double Inv_r_pect[Nregions] = { 0 };
		double max_EMdiff_t = 0; 

		for (int r = 0; r < Nregions; r++) {
			EM_t = EM_t + EM_s[r][t];
			EM0_t = EM0_t + EM0_s[r][t];
			C_Reduce_t = C_Reduce_t + gamma_factor * (CanopyTree[r][t] - CanopyTree[r][0]) * District[r][0] * 2589988 / tree_area;
			Cost_t = Cost_t + Cost_s[r][t];

			EM_r_pect[r] = EM_s[r][t] / (District[r][1] * pow((1 + APG), t)) * 1e5 + 0.0001;   
			EM_pect[r][t] = EM_s[r][t] / (District[r][1] * pow((1 + APG), t)) * 1e5;
			inv_PCap[r][t] = (Cost_Tree_t[r][t] + Cost_CP_t[r][t] + Cost_CR_t[r][t] + Cost_CC_t[r][t]) / (District[r][1] * pow((1 + APG), t));
			Inv_r_pect[r] = (Cost_Tree_t[r][t] + Cost_CP_t[r][t] + Cost_CR_t[r][t] + Cost_CC_t[r][t]) / (District[r][1] * pow((1 + APG), t));
		
			if (EM_r_pect[r] > EM_t_max) {
				EM_t_max = EM_r_pect[r];
			}

			if (EM_r_pect[r] < EM_t_min) {
				EM_t_min = EM_r_pect[r];
			}
			
		}

		// Total EM, Cost, CO2_reduction 
		TEM = TEM + EM_t;
		TEM0 = TEM0 + EM0_t;
		C_Reduce = C_Reduce + C_Reduce_t;
		Cost = Cost + Cost_t;

		// Mortality difference: mean and max
		max_EMdiff_t = EM_t_max - EM_t_min; 
		if ( max_EMdiff_t > max_max_EMdiff) {
			max_max_EMdiff = max_EMdiff_t;
		}     // update the max_EMdiff if the difference is larger than the current max
		
		mean_max_EMdiff = mean_max_EMdiff + max_EMdiff_t / Nyears;  // calculate the mean of one-year max EMdiff

		// Gini: mean and max
		gini_t = GINI_index(EM_r_pect);
		gini_inv = GINI_index(Inv_r_pect);

		if (gini_t > max_gini_t) {
			max_gini_t = gini_t;
		}   //  update the max of gini if the gini is greater than the current max 

		mean_gini_t = mean_gini_t + gini_t / Nyears;    // calculate the mean of one-year gini 

		// One-year mortality: max
		if (EM_t > max_EM) {
			max_EM = EM_t; 
		} // update the max_EM if the difference is larger than the previous one 
		
		// print out results
		fprintf(outputFile, "%d %f %f %f %f %f %f ", (t + 2019), EM0_t, EM_t, Cost_t, gini_t, max_EMdiff_t, C_Reduce_t);  // cost in DPS_invst.out

		for (int r = 0; r < Nregions; r++) {
			fprintf(outputFile, "%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ", \
				HWD0[r][(t + 4)], HWD[r][(t + 4)], EM0_s[r][t], EM_s[r][t], inv_PCap[r][t], \
				Cost_Tree_t[r][t], Cost_CP_t[r][t], Cost_CR_t[r][t], Cost_CC_t[r][t], \
   		        EM0_age64poor[r][t], EM0_age64rich[r][t], EM0_age65poor[r][t], EM0_age65rich[r][t],\
				EM_age64poor[r][t], EM_age64rich[r][t], EM_age65poor[r][t], EM_age65rich[r][t], \
			    pop_age64poor[r][t], pop_age64rich[r][t], pop_age65poor[r][t], pop_age65rich[r][t], \
			    CumTree[r][t], CanopyTree[r][t], CumCP[r][t], CumCR[r][t]);
		}
		fprintf(outputFile, "\n");
	}
	
	REM = TEM0 - TEM;         // annual mortaity reduction 
	Equity = max_max_EMdiff;    // maximum of one-year mortality difference  
	Equity2 = mean_max_EMdiff;
	Equity3 = max_gini_t;        // maximum of one-year gini index 
	Equity4 = mean_gini_t;
	Reliability = max_EM;     // maximum of one-year mortality 

	objs[0] = TEM0;
	objs[1] = TEM;
	objs[2] = REM;            // max
	objs[3] = Cost;
	objs[4] = Equity;
	objs[5] = Equity2;
	objs[6] = Equity3;
	objs[7] = Equity4;
	objs[8] = Reliability;
	objs[9] = C_Reduce;       // max 

	//calculate summer mean temperature 
	for (int i = 1; i < (Nyears + 1); i++) {
		for (int j = 0; j < Nsummer; j++) {
			objs[10] = objs[10] + temp[ClimModel][j][i];
		}
	}
	objs[10] = objs[10] / Nyears / Nsummer;

}

// ===================== Main Function ================================== // 

int main(int argc, char* argv[]) {

	// Main FUN. Part I: INPUT file read

	/*District Mat: 0.Area (sq.mi); 1.Population; 2.Age>=65; 3.Age>=65(%)
	4.Poor; 5.Poor(%); 6.dT; 7.TreePect; 8.RoofPect; 9.RoadPect; 10.TreeAvailPect
	11. ImpSuf NonBuild Pect */

	ifstream infile;
	infile.open("./../../Data/District.txt");
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
	infile1.open("./../../Data/Scenarios.txt");         // use RE_Scenarios.txt for out-of-sample re-evaluation and Set Nsamples = 3300 in line 27
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

	ifstream infile2;
	infile2.open("./../../Data/temp20192039.txt");
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
	infile3.open("./../../Data/BWI_temp.txt");
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
	double vars[nvars];
	double objs[nobjs];

	MOEA_Init(nobjs, 0);

	string soln = argv[1];
	sprintf(output_file_name, "./../DecisionVariables/Output/details/DPSdetail_%s.out", soln.c_str());

	while (MOEA_Next_solution() == MOEA_SUCCESS) {
		MOEA_Read_doubles(nvars, vars);

		outputFile = fopen(output_file_name, "w");
		if (!outputFile) {
			perror("Unable to open final output file\n");
		}

		for (int s = 0; s < Nsamples; s++) {
			//printf("sample %d\n", s);
			ClimModel = (int)ScenMat[s][0]; HWdef = (int)ScenMat[s][1];
			TR_tree = ScenMat[s][2]; TR_CP = ScenMat[s][3]; TR_CR = ScenMat[s][4];
			CCuse = ScenMat[s][5]; p_vage65 = ScenMat[s][6]; p_vpoor = ScenMat[s][7];
			alpha = ScenMat[s][8]; beta = ScenMat[s][9];
			APG = ScenMat[s][10]; AAR = ScenMat[s][11]; APR = ScenMat[s][12]; DisR = ScenMat[s][13];
			HWRR_poor = ScenMat[s][14]; HWRR_age65 = ScenMat[s][15]; HWRR_age64 = ScenMat[s][16];
			CityHEAT(vars, objs, 0);    // unconstrained problem  
			MOEA_Write(objs, NULL);
		}

		fclose(outputFile);

	}

	MOEA_Terminate();

	return EXIT_SUCCESS;
}

