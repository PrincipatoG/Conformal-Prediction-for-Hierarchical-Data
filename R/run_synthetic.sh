#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node 36
#SBATCH --wckey="P1274:R"
#SBATCH --time=3-00:00:00
#SBATCH --job-name="running_simulation"
#SBATCH --error=error.txt
#SBATCH --output=output.txt
#SBATCH --partition=cn
#SBATCH --exclusive

# Activate the environment R
source /projets/datascience_retd/dist/anaconda/bin/activate r-env

echo "Running study.R"
Rscript R/study.R 5 10000 1000 # n, nobs, nsimu
# echo "Running results_analysis.R"
# Rscript R/results_analysis.R