#/bin/bash

#module load gcc python
module load python


python pareto.py ./sets/*.txt -o 110-114 -e 5 1e7 1e-5 10 5e4 --output DPS.resultfile --delimiter=" " --comment="#" 
cut -d ' ' -f 111-115 DPS.resultfile > Overall.reference
