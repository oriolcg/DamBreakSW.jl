using DrWatson
@quickactivate "DamBreakSW"
using TimerOutputs

include(srcdir("DamBreakSW.jl"))

# Warm-up run
params_coarse_2D = DamBreakSW.DamBreak_benchmark_2D_params(
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
params_coarse_2D_ASGS = DamBreakSW.DamBreak_benchmark_2D_params(
  formulation=:ASGS,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
params_coarse_2D_Smagorinsky = DamBreakSW.DamBreak_benchmark_2D_params(
  formulation=:Smagorinsky,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
params_coarse_2D_conservative = DamBreakSW.DamBreak_benchmark_2D_params(
  formulation=:conservative_Galerkin,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
# DamBreakSW.main(params_coarse_2D)
DamBreakSW.main(params_coarse_2D_ASGS)
# DamBreakSW.main(params_coarse_2D_Smagorinsky)
# DamBreakSW.main(params_coarse_2D_conservative)

to = TimerOutput()

# Fine run
params_fine_2D = DamBreakSW.DamBreak_benchmark_2D_params(
  mesh_file=datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
  order=2,
  verbose=false,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_2D",
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=10.0,ρ∞=0.0)
)
params_fine_2D_ASGS = DamBreakSW.DamBreak_benchmark_2D_params(
  mesh_file=datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
  order=2,
  formulation=:ASGS,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_2D_ASGS",
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=10.0,ρ∞=0.0,Δt=0.2)
)
params_fine_2D_Smagorinsky = DamBreakSW.DamBreak_benchmark_2D_params(
  mesh_file=datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
  order=2,
  formulation=:Smagorinsky,
  verbose=false,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_2D_Smagorinsky",
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=10.0,ρ∞=0.0)
)
params_fine_2D_conservative = DamBreakSW.DamBreak_benchmark_2D_params(
  mesh_file=datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
  order=1,
  formulation=:conservative_Galerkin,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_2D_conservative",
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=10.0,ρ∞=0.0)
)
# @timeit to "DamBreak_benchmark_2D" DamBreakSW.main(params_fine_2D)
@timeit to "DamBreak_benchmark_2D_ASGS" DamBreakSW.main(params_fine_2D_ASGS)
# @timeit to "DamBreak_benchmark_2D_Smagorinsky" DamBreakSW.main(params_fine_2D_Smagorinsky)
# @timeit to "DamBreak_benchmark_2D_conservative" DamBreakSW.main(params_fine_2D_conservative)

show(to)
