`cd` into the new  `/Optimization` folder in your HPC. First, run `mkdir sets runtime` to create two folders to save the output from the optimization. 
Ensure that the following files are also in the `Optimization` directory: 
`find_refSets.sh` \\ 
`find_runtime_metrics.sh` \\ 
`get_objs.sh` \\ 
`JobSub.sh`. 
Edit`JobSub.sh` to setup the number of seeds (`line 15`, currently is {1..10} for 10 seeds) and number of function evaluation(`line 17`, currently is 200000 function evaluations) to use. 

Then run `sbatch JobSub.sh` to submit your MOEA job (UrbanHeatMPI from last step) to the HPC and wait for the completion (may take several hours). Once the job is finished, x (x = number of seeds set in the previous step) output files are generated in the `sets` and `runtime` folders with the name of `UrbanHeat_Sx.txt` and `UrbanHeat_Sx.runtime`, where x is the number of seeds. 

Then, you are going to analyze the pareto solutions generated. Download [MOEAFramework-2.13-Demo.jar](https://github.com/MOEAFramework/MOEAFramework/releases/) and [pareto.py](https://github.com/matthewjwoodruff/pareto.py) and upload the two files to this folder. Next, run `mkdir objs`, `sh get_objs.sh`, this will extract the values of objectives and save then in `UrbanHeat_Sx.objs` in the folder objs. Then, run `sh find_refSets.sh`, which will find you a reference set from the current pareto solutions in the `DPS.resultfile` and `Overall.reference`. This reference set will be used to calculate the Hypervolume of generated pareto front. Then, run `mkdir output metrics`, and `sbatch find_runtime_metrics.sh`. `find_runtime_metrics.sh` will generate six indicators of algorithm convergence. It runs multiple jobs in parallel. 

If you want to speed up the process of calculating `run_time_metrics.sh`, you use `WFG2` algorithm (copy the [WFG](https://github.com/MOEAFramework/Hypervolume) folder to your directory). Compile the `WFG2` by using `make` and then copy `wfg2` to your directory (You may need chmod 755 wfg2 to get permission from your machine). Download global.properties from here and copy to your directory. Insert the following two lines in global.properties: 

`org.moeaframework.core.indicator.hypervolume = ./wfg2 {2}` \
`org.moeaframework.core.indicator.hypervolume_inverted = true`

Then just sbatch `find_runtime_metrics.sh`
