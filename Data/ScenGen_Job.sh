#!/bin/bash
#SBATCH -J ScenGenerate
#SBATCH --partition shared 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=rshi8@jhu.edu 

module load R/3.6.1

cd /home-4/rshi8@jhu.edu/scratch/Rui/UrbanHeat/simulation/Data/ 
Rscript --no-save --no-restore ./ScenarioGenerator.R

