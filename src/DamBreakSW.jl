module DamBreakSW

using DrWatson
@quickactivate "DamBreakSW"

using Gridap
using GridapGmsh: gmsh, GmshDiscreteModel
using Parameters

include("time_integrator.jl")
include("weak_form.jl")
include("DamBreak_benchmark_1D.jl")
include("DamBreak_benchmark_2D.jl")
include("DamBreak_building.jl")
# include(raw"C:\Users\ocolomesgene\Progs\gmsh\gmsh-4.13.1-Windows64-sdk\lib\gmsh.jl")
include("mesh_generation.jl")

export DamBreak_benchmark_1D_params
export DamBreak_benchmark_2D_params
export main

end
