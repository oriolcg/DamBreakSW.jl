module Run_Case
using DrWatson
@quickactivate "DamBreakSW"

using MPI
MPI.Init()
comm = MPI.COMM_WORLD
np = MPI.Comm_size(comm)

include(srcdir("DamBreakSW.jl"))

# Run
params_fine_2D = DamBreakSW.DamBreak_building_params(
  mesh_file=datadir("meshes","DamBreak_building_0.0_fine.msh"),
  x₀=0.0,
  ϵ=0.1,
  order=1,
  # order=2,
  verbose=true,
  vtk_output=true,
  Δtout=0.05,
  vtk_folder="DamBreak_building_0.0_fine_ASGS_serial",
  # vtk_folder="DamBreak_building_0.0_fine_Smagorinsky_serial",
  formulation=:ASGS,
  # formulation=:Smagorinsky,
  # ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0,Δt=0.01,T=30.0)
  ode_solver_params=DamBreakSW.ODE_solver_params(:θ_method;Δt=0.01,T=30.0)
  # ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3)
)
# DamBreakSW.main_parallel(np,params_fine_2D)
DamBreakSW.main(params_fine_2D)

MPI.Finalize()
end
