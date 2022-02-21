#!/bin/bash
#SBATCH --job-name=FindMetrics_ma
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p parallel				# Queue name "standard" (serial)
#SBATCH -t 2:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --array=1-10				# Array of jobs to loop through
#SBATCH --mail-type=ALL                  	# email at Begin and End of job

module load java

java -cp MOEAFramework-2.13-Demo.jar org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
	-d 5 -i ./objs/UrbanHeat_S${SLURM_ARRAY_TASK_ID}.obj -r Overall.reference \
	-o ./metrics/UrbanHeat_S${SLURM_ARRAY_TASK_ID}.metrics
