#!/bin/bash
#SBATCH --job-name=MOEA_ma
#SBATCH --nodes=3            					 		        # Total number of nodes to request (up to 120)
#SBATCH --ntasks-per-node=24		               				        # Number of processors per node (up to 20)
#SBATCH --partition=parallel          				 			# Queue name "parallel"
#SBATCH --time=10:00:00      					 			# Run time (hh:mm:ss) - up to 36 hours
#SBATCH --mail-user=rshi8@jhu.edu   		                			# address for email notification
#SBATCH --mail-type=ALL                  	                		        # email at Begin and End of job

#module load gcc openmpi 

#pwd
# Your commands go here
# arguments are <seed> <NFE>
for i in {1..10}
do
  mpirun ../Compile/UrbanHeatMPI $i 200000
done
