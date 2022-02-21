

Result = read.csv("./../../Optimization/DPS.resultfile", sep ="", header = FALSE)

DV = Result[,1:110]
l_DPS = dim(DV)[1]

for (i in 1:l_DPS){
  a = DV[i,]
  write.table(a, file = paste("./DPS", i, ".vars", sep=""), row.names = F, col.names = F)
  
}

