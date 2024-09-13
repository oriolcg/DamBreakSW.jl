module DamBreakSW

using DrWatson
@quickactivate "DamBreakSW"

using Gridap
using GridapGmsh: gmsh, GmshDiscreteModel
using GridapDistributed
using GridapPETSc
using PartitionedArrays
using Parameters

include("time_integrator.jl")
include("weak_form.jl")
include("DamBreak_benchmark_1D.jl")
include("DamBreak_benchmark_2D.jl")
include("DamBreak_building.jl")
# include(raw"C:\Users\ocolomesgene\Progs\gmsh\gmsh-4.13.1-Windows64-sdk\lib\gmsh.jl")
include("mesh_generation.jl")

export DamBreak_benchmark_1D_params
export DamBreak_benchmark_2D_params
export main, main_parallel

options_mumps = "-snes_type newtonls \
-snes_linesearch_type basic  \
-snes_linesearch_damping 1.0 \
-snes_rtol 1.0e-6 \
-snes_atol 1.0e-8 \
-snes_max_it 20 \
-ksp_error_if_not_converged true \
-ksp_converged_reason -ksp_type preonly \
-pc_type lu \
-pc_factor_mat_solver_type mumps \
-mat_mumps_icntl_7 0 \
-mat_mumps_icntl_14 500000"

function main_parallel(np,params)
  current_path = pwd()
  cd(datadir("sims",params.vtk_folder))
  with_mpi() do distribute
    options = options_mumps
    ranks = distribute(LinearIndices((np,)))
    GridapPETSc.with(args=split(options)) do
      main(ranks,params)
    end
  end
  cd(current_path)
end

end
