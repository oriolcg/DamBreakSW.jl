#!/bin/bash

#SBATCH --job-name="compile_DamBreakSW"
#SBATCH -p thin
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -o stdout
#SBATCH -e stderr

source modules_snellius.sh
srun julia --project=../ -O0 --check-bounds=no --color=yes compile.jl

