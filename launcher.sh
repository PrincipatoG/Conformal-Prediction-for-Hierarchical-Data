#!/bin/bash

#SBATCH -N 1
#SBATCH --wckey="P1274:R"
#SBATCH --time=3-00:00:00
#SBATCH --job-name="launching_simulation"
#SBATCH --output=launcher_output.txt
#SBATCH --error=launcher_error.txt
#SBATCH --partition=cn
#SBATCH --enable-turbo

# Activate the environment R
source /projets/datascience_retd/dist/anaconda/bin/activate r-env

# configs=(
#   "3 12 1000000 1000 1"
#   "4 12 1000000 1000 2"
#   "3 144 1000000 1000 3"
#   "4 144 1000000 1000 4"
#   "3 1728 1000000 1000 5"
#   "4 1728 1000000 1000 6"
# )
configs=(
   "3 12 1000000 1000 21"
   "4 12 1000000 1000 22"
)

for i in "${!configs[@]}"; do
  config_number=$((i+1))
  echo "${configs[$i]}" 
  sbatch --wait -p cn 'run.sh' "${configs[$i]}"
  srun --exclusive --ntasks=1 Rscript R/gather.R --config "${configs[$i]}"
done

srun --exclusive --ntasks=1 Rscript R/LaTeX_Table.R
