#!/bin/bash

#SBATCH --wckey=xxxxxxx
#SBATCH --time=3-00:00:00
#SBATCH --job-name="running_simulation"
#SBATCH --array=0-999  
#SBATCH --output=logs/run_output_%A_%a.txt
#SBATCH --error=logs/run_error_%A_%a.txt

# Activate the environment R
source xxxxxxx r-env

CONFIG_FILE=$1
srun Rscript R/study.R --simulation_index $SLURM_ARRAY_TASK_ID --config "$CONFIG_FILE"
