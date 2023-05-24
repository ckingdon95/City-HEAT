#!/bin/bash

#SBATCH --job-name=MOEA_ma
#SBATCH --account=fc_anthofflab
#SBATCH --partition=savio          				 			
#SBATCH --nodes=1            					 		        
#SBATCH --cpus-per-task=24		               				        
#SBATCH --time=1:00:00      					 			

#module load gcc openmpi 

#pwd
# Your commands go here
# arguments are <seed> <NFE>
for i in {1..10}
do
  mpirun ../Compile/UrbanHeatMPI $i 200000
done
