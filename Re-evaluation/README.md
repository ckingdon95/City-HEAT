1. `cd` to `/Re-evaluation` folder in your HPC. Run `mkdir DecisionVariables` and `cd` to this folder. 

3. Run `Rscript ResultExtract.R` to generate multiple `.vars` files, each of which contains the values of decision variables for a specific solution.

5. Run `mkdir Output` and `cd` to `Output`, run `mkdir details`.

7. `cd` back to `/Re-evaluation`.


9. Run `mkdir DetailOut` and `cd` to this folder.


11. Download `moeaframework.cpp` and `moeaframework.h` [here](http://moeaframework.org/). 


13. Upload `ReEvlDetail_Array.sh`, `makefile`, `UHMOEA_ReEvl.cpp` to this folder.


15. Run `make` and `UrbanHeat_RE_detail.exe` will be generated. 


17. Next, `sbatch ReEvlDetail_Array.sh`. This will generate detailed Re-evaluation results of the pareto solutions of your interests (in “/Re-evaluation/DecisionVariables” folder).                        
    + The Re-Evaluation results are saved in `/Re-evaluation/DecisionVariables/Output` (re-simulated objective functions) and “/Re-evaluation/DecisionVariables/Output/details” (other interesting metrics, more detailed resolutions, e.g., investment in each year), respectively.


18. `cd` back to `/Re-evaluation` and upload `Robustness.R`.


20. Run `Rscript Robustness.R` to calculate robustness metrics for each solution. The result is generated and saved in `Robustness.out`.
