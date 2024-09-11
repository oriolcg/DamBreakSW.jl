module DamBreakSW

using DrWatson
@quickactivate "DamBreakSW"

using Gridap
using GridapGmsh
using Parameters

include("time_integrator.jl")
include("weak_form.jl")
include("DamBreak_benchmark_1D.jl")
include("DamBreak_benchmark_2D.jl")
include("DamBreak_building.jl")

export DamBreak_benchmark_1D_params
export DamBreak_benchmark_2D_params
export main

end
