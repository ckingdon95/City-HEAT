#!/bin/bash
#SBATCH --job-name=MOEA_ma
#SBATCH --nodes=1            					 		        # Total number of nodes to request (up to 120)
#SBATCH --cpus-per-task=24		               				        # Number of processors per node (up to 20)
#SBATCH --partition=savio2          				 			# Queue name "parallel"
#SBATCH --time=10:00:00      					 			# Run time (hh:mm:ss) - up to 36 hours
#SBATCH --account=fc_anthofflab

#module load gcc openmpi 

#pwd
# Your commands go here
# arguments are <seed> <NFE>
for i in {1..10}
do
  mpirun ../Compile/UrbanHeatMPI $i 200000
done
