`Data` folder contains all inputs used in the paper for Baltimore. Users need to prepare these datasets for the city of thier interests.

1. `District.txt` contains information for 11 Planning Districts in Baltimore City, each column means:

  * column 1.  Area (sq.mi) of district 
  * column 2.  Total population 
  * column 3.  Number of people over 65 years old  
  * column 4.  Percentage of people over 65 years old (%)
  * column 5.  Number of people below poverty line
  * column 6.  Percentage of people below poverty line (%)
  * column 7.  Temperature deviation of district 
  * column 8.  Tree canopy coverage (%) 
  * column 9.  Roof area percentage (%) 
  * column 10.  Rood area percentage (%)
  * column 11. Grass lands for tree planting (%)
  * column 12. Impervious surface (excluding Road & Roof) area percentage (%)

Demographic data from [American Community Survey]( https://www.census.gov/programs-surveys/acs) (columns 2-6) \
Temperature data from [JHU monitoring practice]( https://journals.ametsoc.org/view/journals/apme/56/1/jamc-d-16-0232.1.xml) (column 7) \
Land use data from [Chesapeake Phase 6 land use data]( https://www.chesapeakeconservancy.org/conservation-innovation-center/high-resolution-data/land-use-data-project/) (columns 1 & 8-12) \

Each row is for disctrt `Central`, `Downtown`, `East`, `North`, `Northeast`, `Northwest`, `South`, `Southeast`, `Southwest`, `Southwest Partnership`, and `West` respectively. 


2. `Scenarios.txt` and `RE_Scenarios.txt` contain the values of uncertain parameters for 1,500 optimization scenarios and 3,300 re-evaluation scenarios. 
  * column 1. Climate model	id
  * column 2. Heat wave definition id	
  * column 3. Temperature reduction effect for tree
  * column 4. Temperature reduction effect for cool pavement	
  * column 5. Temperature reduction effect for cool roof	
  * column 6. Cooling center usage rate	
  * column 7. Cooling center vistors over 65 years old (%)	
  * column 8. Cooling center vistors below poverty line (%)	
  * column 9. Percentage of impervious surface being converted to tree (%)
  * column 10. Spillover effect from adjecent districts	
  * column 11. Average annual population grouwth rate (%)	
  * column 12. Average annual aging rate (%)	
  * column 13. Average annual poverty rate (%)	
  * column 14. Real discount rate	
  * column 15. Heat wave relative mortality risk for people below poverty line	
  * column 16. Heat wave relative mortality risk for people over 65 years old
  * column 17. Heat wave relative mortality risk for people below 65 years old

3. `temp20202039.txt` contains 32 20-year daily temperature projections (32*20*365 data points in total) retrieved from [NA-CORDEX dataset](https://na-cordex.org/) and [LOCA dataset](http://loca.ucsd.edu/). 

	GCM	RCM	RCP	Bias-correction	Source
1	CanESM2	CanRCM4	rcp85	mbcn-Daymet	NA-Cordex
2	CanESM2	CRCM5-UQAM	rcp85	mbcn-Daymet	NA-Cordex
3	GEMatm-Can	CRCM5-UQAM	rcp85	mbcn-Daymet	NA-Cordex
4	GEMatm-MPI	CRCM5-UQAM	rcp85	mbcn-Daymet	NA-Cordex
5	GFDL-ESM2M	RegCM4	rcp85	mbcn-Daymet	NA-Cordex
6	GFDL-ESM2M	WRF	rcp85	mbcn-Daymet	NA-Cordex
7	HadGEM2-ES	WRF	rcp85	mbcn-Daymet	NA-Cordex
8	MPI-ESM-LR	CRCM5-UQAM	rcp85	mbcn-Daymet	NA-Cordex
9	MPI-ESM-LR	RegCM4	rcp85	mbcn-Daymet	NA-Cordex
10	MPI-ESM-LR	WRF	rcp85	mbcn-Daymet	NA-Cordex
11	MPI-ESM-MR	CRCM5-UQAM	rcp85	mbcn-Daymet	NA-Cordex
12	GFDL-ESM2M	WRF	rcp85	mbcn-gridMET	NA-Cordex
13	HadGEM2-ES	WRF	rcp85	mbcn-gridMET	NA-Cordex
14	MPI-ESM-LR	CRCM5-UQAM	rcp85	mbcn-gridMET	NA-Cordex
15	MPI-ESM-LR	RegCM4	rcp85	mbcn-gridMET	NA-Cordex
16	MPI-ESM-LR	WRF	rcp85	mbcn-gridMET	NA-Cordex
17	MPI-ESM-MR	CRCM5-UQAM	rcp85	mbcn-gridMET	NA-Cordex
18	CanESM2	CanRCM4	rcp85	mbcn-gridMET	NA-Cordex
19	CanESM2	CRCM5-UQAM	rcp85	mbcn-gridMET	NA-Cordex
20	GEMatm-Can	CRCM5-UQAM	rcp85	mbcn-gridMET	NA-Cordex
21	GEMatm-MPI	CRCM5-UQAM	rcp85	mbcn-gridMET	NA-Cordex
22	GFDL-ESM2M	RegCM4	rcp85	mbcn-gridMET	NA-Cordex
23	CanESM2	Statistical	rcp45	none	LOCA
24	GFDL-ESM2M	Statistical	rcp45	none	LOCA
25	HadGEM2-ES	Statistical	rcp45	none	LOCA
26	MPI-ESM-LR	Statistical	rcp45	none	LOCA
27	MPI-ESM-MR	Statistical	rcp45	none	LOCA
28	CanESM2	Statistical	rcp85	none	LOCA
29	GFDL-ESM2M	Statistical	rcp85	none	LOCA
30	HadGEM2-ES	Statistical	rcp85	none	LOCA
31	MPI-ESM-LR	Statistical	rcp85	none	LOCA
32	MPI-ESM-MR	Statistical	rcp85	none	LOCA

4. `BWI_temp.txt` contains temperature observations at BWI airport (station: USW00093721)
). 
