using DrWatson
@quickactivate "DamBreakSW"
using TimerOutputs

include(srcdir("DamBreakSW.jl"))

# Warm-up run
params_coarse_1D = DamBreakSW.DamBreak_benchmark_1D_params(
  order=2,
  verbose=false,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
params_coarse_1D_ASGS = DamBreakSW.DamBreak_benchmark_1D_params(
  order=2,
  formulation=:ASGS,
  verbose=false,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
params_coarse_1D_Smagorinsky = DamBreakSW.DamBreak_benchmark_1D_params(
  order=2,
  formulation=:Smagorinsky,
  verbose=false,
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;ρ∞=0.0)
)
params_coarse_1D_Smagorinsky_EXRK = DamBreakSW.DamBreak_benchmark_1D_params(
  order=2,
  formulation=:Smagorinsky,
  verbose=false,
  ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3)
)
# DamBreakSW.main(params_coarse_1D)
DamBreakSW.main(params_coarse_1D_ASGS)
# DamBreakSW.main(params_coarse_1D_Smagorinsky)
# DamBreakSW.main(params_coarse_1D_Smagorinsky_EXRK)

to = TimerOutput()

# Fine run
params_1D = DamBreakSW.DamBreak_benchmark_1D_params(
  nₓ=400,
  order=2,
  verbose=false,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_1D",
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=60.0,ρ∞=0.0)
)
params_1D_ASGS = DamBreakSW.DamBreak_benchmark_1D_params(
  nₓ=400,
  order=2,
  formulation=:ASGS,
  verbose=false,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_1D_ASGS",
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=60.0,ρ∞=0.0)
)
params_1D_Smagorinsky = DamBreakSW.DamBreak_benchmark_1D_params(
  nₓ=400,
  order=2,
  formulation=:Smagorinsky,
  verbose=false,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_1D_Smagorinsky",
  ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α;T=60.0,ρ∞=0.0)
)
params_1D_Smagorinsky_EXRK = DamBreakSW.DamBreak_benchmark_1D_params(
  nₓ=400,
  order=2,
  formulation=:Smagorinsky,
  verbose=false,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_1D_Smagorinsky_EXRK",
  ode_solver_params=DamBreakSW.ODE_solver_params(:EXRK_SSP_3_3;T=60.0,Δt=0.2)
)
# @timeit to "DamBreak_benchmark_1D" DamBreakSW.main(params_1D)
@timeit to "DamBreak_benchmark_1D_ASGS" DamBreakSW.main(params_1D_ASGS)
# @timeit to "DamBreak_benchmark_1D_Smagorinsky" DamBreakSW.main(params_1D_Smagorinsky)
# @timeit to "DamBreak_benchmark_1D_Smagorinsky_EXRK" DamBreakSW.main(params_1D_Smagorinsky_EXRK)

show(to)
