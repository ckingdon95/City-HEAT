
load("LOCACordex_Proj20192039.RData")
#load("D:/OneDrive - Johns Hopkins University/ResearchProject/UrbanHeat/Year3/RDM/MOEA/data/LOCACordex_Proj20192039.RData")

#temp20192039_1 = round(temp20192039, 4)
#write.table(temp20192039_1, file = "D:/OneDrive - Johns Hopkins University/ResearchProject/UrbanHeat/Year3/RDM/MOEA/Code/temp20192039.txt", col.names = F, row.names = F)


library(lhs)

SampleGenerator <- function(ClimMod.id, n.sample, n.parm, SEED){
  
  climatemodel = data.frame(ClimMod.id)   # 0: climate model
  HWDef = data.frame(1:3)                 # 1: heat wave definition 
  set.seed(SEED)
  LHC.sample = randomLHS(n.sample, n.parm)
  LHSpara = array(0, c(n.sample, n.parm))
  
  #LHS parameter
  df.range = rbind(
    TR.tree.range = c(1.5, 7),                     # 2.  temperature reduction per 1% tree canopy increase 
    TR.CP.range = c(0.25, 1.5),                    # 3.  temperature reduction per 1% cool pavement increase
    TR.CR.range = c(0.2, 1.2),                     # 4.  temperature reduction per 1% cool roof increase  
    cc.UR.range = c(0.2, 0.8)  ,                   # 5.  cooling center utilization rate 
    p.cc.age65.range = c(0.2, 0.8),                # 6.  % of visitors older than 65
    p.cc.poor.range = c(0.5, 0.9)  ,               # 7.  % of visitors below poverty line
    alpha.range = c(0.05, 0.3)      ,              # 8.  % of impervious non-building surfaces can be converted to tree
    beta.range = c(0, 0.5)           ,             # 9.  % of spillover effect to adjacent districts 
    APG.range = c(-1.8/100, 0.62/100) ,            # 10. annual population growth rate
    AAR.range = c(0.1/100, 0.6/100)    ,           # 11. annual population aging rate 
    APR.range = c(0, 0.5/100)           ,          # 12. annual poverty rate  
    DisR.range = c(0.01,0.08)            ,         # 13. discount rate 
    HWRR.poor.range = c(1, 1.1)           )        # 14. odds that poor die
  
  for (i in 1:n.parm){
    LHSpara[,i] = qunif(LHC.sample[,i], df.range[i, 1], df.range[i, 2])
  }
  
  ScenarioMat = merge(climatemodel, HWDef)
  ScenarioMat = merge(ScenarioMat, LHSpara)
  
  hwrr64 = hwrr65 = array(0, dim(ScenarioMat)[1])
  hw_id1 = (ScenarioMat[,2]==1)
  LHC.sample1 = randomLHS(sum(hw_id1==1), 1)
  hwrr64[hw_id1] = hwrr65[hw_id1] = qunif(LHC.sample1, 0.06, 0.2)
  
  hw_id2 = (ScenarioMat[,2]==2)
  LHC.sample2 = randomLHS(sum(hw_id2==1), 1)
  hwrr64[hw_id2] = hwrr65[hw_id2] = qunif(LHC.sample2, 0.84/100, 2.85/100)
  
  hw_id3 = (ScenarioMat[,2]==3)
  LHC.sample3 = randomLHS(sum(hw_id3==1), 1)
  hwrr64[hw_id3] = qunif(LHC.sample3, 5.9/100, 11.2/100)
  hwrr65[hw_id3] = qunif(LHC.sample3, 7.8/100, 14.2/100)
  
  ScenarioMat = cbind(ScenarioMat,hwrr65,hwrr64)
  colnames(ScenarioMat) <- c(
    "ClimateModel", "HWDef", "TR_Tree", "TR_CP", "TR_CR",	"CC_usage",
    "P_vage65",	"P_vpoor", "alpha", "beta",	"APG",	"AAR", "APR",	"DistR",	
    "HWRR_poor",	"HWRR65",	"HWRR64"
  )
  
  return(ScenarioMat)
      
}
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

set.seed(1)
climate.id = sample(1:length(unique_model1), 10) 
climate.id.test = (1:length(unique_model1))[-climate.id] -1
climate.id = climate.id -1
## training samples
MOEA_scenarios = SampleGenerator(climate.id, 50, 13, 1)
MOEA_scenarios = round_df(MOEA_scenarios, 4)
write.table(MOEA_scenarios, file = "Scenarios.txt", row.names = F, col.names = F)
## testing samples
ReEvl_scenarios = SampleGenerator(climate.id.test, 50, 13, 2)
ReEvl_scenarios = round_df(ReEvl_scenarios, 4)
write.table(ReEvl_scenarios, file = "RE_Scenarios.txt", row.names = F, col.names = F)


