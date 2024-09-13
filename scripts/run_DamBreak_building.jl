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
  vtk_folder="DamBreak_building_0.0_coarse_Smagorinsky",
  formulation=:Smagorinsky,
  # ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
  ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3)
)
@timeit to "Coarse" DamBreakSW.main(params_coarse_2D)

params_fine_2D = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_fine.msh"),
  x₀=0.0,
  order=2,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_building_0.0_fine_Smagorinsky",
  formulation=:Smagorinsky,
  # ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0,T=30.0,Δt=0.01)
  ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3;T=30.0,Δt=0.001)
)
@timeit to "Fine" DamBreakSW.main(params_fine_2D)

show(to)
