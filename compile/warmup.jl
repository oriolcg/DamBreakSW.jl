using DrWatson
@quickactivate "DamBreakSW"
include(srcdir("DamBreakSW.jl"))
# Warm-up run
params_coarse_2D = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_coarse.msh"),
  x₀=0.0,
  order=2,
  verbose=true,
  vtk_output=true,
  vtk_folder="DamBreak_building_0.0_coarse_Smagorinsky",
  formulation=:Smagorinsky,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
  # ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3)
)
DamBreakSW.main_parallel(1,params_coarse_2D)
