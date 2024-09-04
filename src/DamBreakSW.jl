module DamBreakSW

using DrWatson
@quickactivate "DamBreakSW"

using Gridap
using GridapGmsh
using Parameters

include("DamBreak_benchmark_1D.jl")

export DamBreak_benchmark_1D_params
export main

end
