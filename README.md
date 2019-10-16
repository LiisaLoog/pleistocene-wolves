This folder contains code and scripts to preform analyses described in Loog et al. (2019) “Modern wolves trace their origin to a late Pleistocene expansion from Beringia.” Molecular Ecology


The 4 directories numbered 1 to 4 correspond to the 4 different types of scenarios described in the manuscript. Each directory contains the C++ code to perform simulations under a given scenario and output files containing simulated summary statistics for ABC based on TMRC’s between demes, calculated as described in the Materials & Methods section of Loog et. al. (2019). The code also performs filter for ABC based on a cut-off of R2 > 90% between simulated and empirical summary statistics. To compile each simulation, please use the g++ command line given as a comment in the top of each simulation code file.

The Data directory contains data files that are used by the C++ simulation code (common to all scenarios) and empirical data used for ABC analysis.

To recreate the summary statistics for ABC for all scenarios, please first compile the simulations programmes in each directory as described above. Then run the bash script do_runs.sh (this could take several days).

To calculate the posterior probabilities and Bayes factors using ABC from simulated summary statistics, please follow these steps: 

1. Run the R scripts for each scenario to generate input files for ABCestimator (files starting with gen.ABC.input.code, e.g. gen.ABC.input.code.1.static.r), making sure the working directory of R is set to the same directory as this readme file. 

2. Install ABCtoolbox (e.g. from github, https://github.com/sakeel/ABCtoolbox.git) and place a copy of the ABCestimator binary in the same directory as this readme file. 

3. Run the bash script run_abc.sh. This will run ABCestimator on all the input files (including the params.txt file containing the parameter settings for the ABC estimator) and extract information for calculating posterior likelihoods. This information is stored in text files in each scenario directory with file names beginning with ABC_RESULTS.

4. Update the table in the file ABC_Marginal_Densities.csv with the new marginal densities and line counts in the ABC_RESULTS files (this step is currently manual). Then use the R script ABC_BF_script.r to calculate posterior probabilities and Bayes factors of each scenario, which were reported in the main text of the manuscript and used to generate figure 1.

Finally, there is one R script for each scenario with file names starting with plot.ABC.output (e.g. plot.ABC.output.code.1.static.r). Use the scripts to generate plots of the posterior parameter distributions for all parameters  (as in supplementary figures 8 and 9).
