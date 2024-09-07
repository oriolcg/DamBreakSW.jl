using DrWatson, Test
@quickactivate "DamBreakSW.jl"
using TimerOutputs
using DataStructures

include(srcdir("DamBreakSW.jl"))

# Set up test parameters
param_dict = OrderedDict(
  :coarse_1D => DamBreakSW.DamBreak_benchmark_1D_params(),
  :coarse_1D_Gα => DamBreakSW.DamBreak_benchmark_1D_params(
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :coarse_1D_Gα_ASGS => DamBreakSW.DamBreak_benchmark_1D_params(
    formulation=:ASGS,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :coarse_1D_Gα_Smagorinsky => DamBreakSW.DamBreak_benchmark_1D_params(
    formulation=:Smagorinsky,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :fine_1D => DamBreakSW.DamBreak_benchmark_1D_params(
    nₓ=400,
    verbose=true
  ),
  :fine_1D_Gα => DamBreakSW.DamBreak_benchmark_1D_params(
    nₓ=400,
    verbose=true,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :fine_1D_Gα_ASGS => DamBreakSW.DamBreak_benchmark_1D_params(
    nₓ=400,
    verbose=true,
    formulation=:ASGS,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :fine_1D_Gα_Smagorinsky => DamBreakSW.DamBreak_benchmark_1D_params(
    nₓ=400,
    verbose=true,
    formulation=:Smagorinsky,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :coarse_2D => DamBreakSW.DamBreak_benchmark_2D_params(),
  :fine_2D => DamBreakSW.DamBreak_benchmark_2D_params(
    mesh_file = datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
    verbose=true
  ),
  :coarse_2D_Gα_ASGS => DamBreakSW.DamBreak_benchmark_2D_params(
    formulation=:ASGS,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :coarse_2D_Gα_Smagorinsky => DamBreakSW.DamBreak_benchmark_2D_params(
    formulation=:Smagorinsky,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  ),
  :fine_2D_Gα_ASGS => DamBreakSW.DamBreak_benchmark_2D_params(
    mesh_file = datadir("meshes","DamBreak_benchmark_2D_fine.msh"),
    verbose=true,
    formulation=:ASGS,
    ode_solver_params=DamBreakSW.ODE_solver_params(:Generalized_α)
  )
)

# Run test suite
println("Starting tests")
to = TimerOutput()
@testset "DamBreakSW.jl tests" begin
  for (key,value) in param_dict
    @timeit to "$key test" begin
      @test DamBreakSW.main(value) == nothing
    end
  end
end

show(to; compact=true)
