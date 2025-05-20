#!/bin/bash

#SBATCH -N 1
#SBATCH --wckey= xxxxxxx
#SBATCH --time=3-00:00:00
#SBATCH --job-name="eval_simulation"
#SBATCH --output=eval_output.txt
#SBATCH --error=eval_error.txt
#SBATCH --partition=cn
#SBATCH --enable-turbo

# Activate the environment R
source xxxxxxx r-env

configs=(
  "3 12 1000000 1000 1"
  "4 12 1000000 1000 2"
  "3 144 1000000 1000 3"
  "4 144 1000000 1000 4"
  "3 1728 1000000 1000 5"
  "4 1728 1000000 1000 6"
)

for i in "${!configs[@]}"; do
  config_number=$((i+1))
  echo "${configs[$i]}" 
  srun --ntasks=1 Rscript R/gather.R --config "${configs[$i]}"
  srun --ntasks=1 Rscript R/plot.R --config "${configs[$i]}"
done

srun --ntasks=1 Rscript R/LaTeX_Table.R
