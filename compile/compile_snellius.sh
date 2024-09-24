#!/bin/bash

#SBATCH --job-name="compile_DamBreakSW"
#SBATCH -p rome
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -o stdout
#SBATCH -e stderr

source modules_snellius.sh
julia --project=../ -e 'using Pkg; Pkg.build("MPI"); Pkg.build("GridapPETSc")'
srun julia --project=../ -O0 --check-bounds=no --color=yes compile.jl

