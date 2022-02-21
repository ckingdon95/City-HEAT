#!/bin/bash
#SBATCH --job-name=RE_ReEvl_ma
#SBATCH --ntasks=1				           # Number of tasks per serial job (must be 1)
#SBATCH --partition=lrgmem				           # Queue name "standard" (serial)
#SBATCH -t 10:00:00				           # Run time per serial job (hh:mm:ss)
#SBATCH --array=1-166			                   # Array of jobs to loop through
#SBATCH --mail-user=rshi8@jhu.edu   		           # address for email notification
#SBATCH --mail-type=ALL                    		   # email at Begin and End of job



./UrbanHeat_RE_detail.exe ${SLURM_ARRAY_TASK_ID} < ../DecisionVariables/DPS${SLURM_ARRAY_TASK_ID}.vars \
						> ../DecisionVariables/Output/DPS${SLURM_ARRAY_TASK_ID}.out

