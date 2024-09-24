#!/bin/sh
#
#SBATCH --job-name="DamBreakSW"
#SBATCH --partition=rome
#SBATCH --time=01:00:00
#SBATCH -n 32
#SBATCH -o stdout-benchmark/slurm-%j-%4t-%n.out
#SBATCH -e stdout-benchmark/slurm-%j-%4t-%n.err

source ../compile/modules_snellius.sh
# julia --project=../ -e 'using Pkg; Pkg.build("MPI"); Pkg.build("GridapPETSc")'
srun julia --project=../ -J ../DamBreakSW.so -O0 --check-bounds=no -e 'include("run_case.jl")'
