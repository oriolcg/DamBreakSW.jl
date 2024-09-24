using DrWatson
@quickactivate "DamBreakSW"
using TimerOutputs

include(srcdir("DamBreakSW.jl"))

to = TimerOutput()

# Warm-up run
params_coarse_2D = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_coarse.msh"),
  x₀=0.0,
  order=1,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_building_0.0_coarse_ASGS",
  formulation=:ASGS,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
  # ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3)
)
params_coarse_2D_ASGS = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_coarse.msh"),
  x₀=0.0,
  order=2,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_building_0.0_coarse_ASGS",
  formulation=:ASGS,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
  # ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3)
)
# @timeit to "Coarse" DamBreakSW.main(params_coarse_2D)
@timeit to "Coarse" DamBreakSW.main(params_coarse_2D_ASGS)

params_medium_2D = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_medium.msh"),
  x₀=0.0,
  order=1,
  verbose=true,
  vtk_output=true,
  Δtout = 0.0,
  vtk_folder="DamBreak_building_0.0_medium_conservative",
  formulation=:conservative_Galerkin,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0,T=0.05,Δt=0.01)
  # ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3;T=30.0,Δt=0.001)
)
params_medium_2D_ASGS = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_medium.msh"),
  x₀=0.0,
  order=2,
  verbose=true,
  vtk_output=true,
  Δtout = 0.0,
  vtk_folder="DamBreak_building_0.0_medium_ASGS",
  formulation=:ASGS,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0,T=1.0,Δt=0.01)
  # ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3;T=30.0,Δt=0.001)
)
# @timeit to "Medium" DamBreakSW.main(params_medium_2D)
@timeit to "Medium" DamBreakSW.main(params_medium_2D_ASGS)

show(to)
