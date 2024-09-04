using DrWatson
@quickactivate "DamBreakSW"

include(srcdir("DamBreakSW.jl"))

# Warm-up run
params_1D = DamBreakSW.DamBreak_benchmark_1D_params()
DamBreakSW.main(params_1D)

# Fine run
params_1D = DamBreakSW.DamBreak_benchmark_1D_params(
  nₓ=400,
  order=2,
  T = 60.0,
  verbose=false,
  vtk_output=true,
  vtk_folder="DamBreak_benchmark_1D"
)
DamBreakSW.main(params_1D)
