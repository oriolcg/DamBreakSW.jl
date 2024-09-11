using DrWatson
@quickactivate "DamBreakSW"
using TimerOutputs

include(srcdir("DamBreakSW.jl"))

to = TimerOutput()

# Warm-up run
params_coarse_2D = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_coarse.msh"),
  x₀=0.0,
  order=2,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_building_0.0_coarse_ASGS",
  formulation=:ASGS,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
@timeit to "Coarse" DamBreakSW.main(params_coarse_2D)

params_fine_2D = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_fine.msh"),
  x₀=0.0,
  order=2,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_building_0.0_fine_ASGS",
  formulation=:ASGS,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0,T=15.0,Δt=0.01)
)
@timeit to "Fine" DamBreakSW.main(params_fine_2D)

# # Fine run
# params_fine_2D = DamBreakSW.DamBreak_benchmark_2D_params(
#   mesh_file=datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
#   order=2,
#   verbose=false,
#   vtk_output=true,
#   vtk_folder="DamBreak_benchmark_2D",
#   ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=10.0,ρ∞=0.0)
# )
# params_fine_2D_ASGS = DamBreakSW.DamBreak_benchmark_2D_params(
#   mesh_file=datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
#   order=2,
#   formulation=:ASGS,
#   verbose=false,
#   vtk_output=true,
#   vtk_folder="DamBreak_benchmark_2D_ASGS",
#   ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=10.0,ρ∞=0.0)
# )
# params_fine_2D_Smagorinsky = DamBreakSW.DamBreak_benchmark_2D_params(
#   mesh_file=datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
#   order=2,
#   formulation=:Smagorinsky,
#   verbose=false,
#   vtk_output=true,
#   vtk_folder="DamBreak_benchmark_2D_Smagorinsky",
#   ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=10.0,ρ∞=0.0)
# )
# @timeit to "DamBreak_benchmark_2D" DamBreakSW.main(params_fine_2D)
# # @timeit to "DamBreak_benchmark_2D_ASGS" DamBreakSW.main(params_fine_2D_ASGS)
# @timeit to "DamBreak_benchmark_2D_Smagorinsky" DamBreakSW.main(params_fine_2D_Smagorinsky)

show(to)
