# City-HEAT

City-Heat Equity Adaptation Tool (C++, R, Bash)

This repository contains all of the code used for the study described in the following working paper:


Shi, R., Hobbs, B., Quinn, D., Knopman, D., Lempert, R. (2021). City-Heat Equity Adaptation Tool (City-HEAT): A Comprehensive Decision Support Tool Considering Multiple Objectives, Uncertainty, and Adaptive Management for Urban Heat Adaptation

Intended for use with the MOEAFramework, Borg and pareto.py. Licensed under the GNU Lesser General Public License.

This research project was conducted using computational resources (and/or scientific computing services) at the Maryland Advanced Research Computing Center ([MARCC](https://www.marcc.jhu.edu/)). The MARCC high-performance computing (HPC) environment utilizes the [SLURM](https://www.marcc.jhu.edu/getting-started/running-jobs/) queuing system for its parallel programming tasks. This code has _not_ been tested nor run on a personal machine.

Contents:

`Data`: Directory containing all of the data input needed in this paper, including district-level demographic and temperature deviation data, temperature projections from 2020 to 2039, temperature observations at BWI airport, and parameter values in each scenario for optimization and re-evaluation, respectively.

`Compile`: Directory containing code for Optimization. 

`Optimization`: Directory containing the code for optimizing the City-HEAT using DPS and Borg-MOEA

`Re-evaluation`: Directory containing the code for re-evaluating the policies found from the above optimization on a set of alternative SOWs and calculating the satisficing metric for each policy

To reproduce the results from this study, follow the steps given in the `README` of the `Data`, `Compile`, `Optimization`, and `Re-evaluation` directory.
