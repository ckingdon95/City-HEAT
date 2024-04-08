#!/bin/bash

#SBATCH --job-name=RE_ReEvl_ma
#SBATCH --account=fc_anthofflab
#SBATCH --partition=savio_normal          				 			
#SBATCH --nodes=1            					 		        
#SBATCH --cpus-per-task=24		               				        
#SBATCH --time=1:00:00   
#SBATCH --array=1-166			                   # Array of jobs to loop through
#SBATCH --mail-user=ckingdon@berkeley.edu  		           # address for email notification
#SBATCH --mail-type=ALL                    		   # email at Begin and End of job



./UrbanHeat_RE_detail.exe ${SLURM_ARRAY_TASK_ID} < ../DecisionVariables/DPS${SLURM_ARRAY_TASK_ID}.vars \
						> ../DecisionVariables/Output/DPS${SLURM_ARRAY_TASK_ID}.out

