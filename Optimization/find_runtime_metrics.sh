#!/bin/bash
#SBATCH --job-name=FindMetrics_ma
#SBATCH --account=fc_anthofflab
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p parallel				# Queue name "standard" (serial)
#SBATCH --array=1-3				# Array of jobs to loop through
#SBATCH --partition=savio          				 			
#SBATCH --nodes=1            					 		        
#SBATCH --time=1:00:00      

module load java

java -cp MOEAFramework-2.13-Demo.jar org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
	-d 5 -i ./objs/UrbanHeat_S${SLURM_ARRAY_TASK_ID}.obj -r Overall.reference \
	-o ./metrics/UrbanHeat_S${SLURM_ARRAY_TASK_ID}.metrics
