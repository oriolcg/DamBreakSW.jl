using DrWatson, Test
@quickactivate "DamBreakSW.jl"
using TimerOutputs

include(srcdir("DamBreakSW.jl"))

# Run test suite
println("Starting tests")
to = TimerOutput()

@testset "DamBreakSW.jl tests" begin
  params_1D = DamBreakSW.DamBreak_benchmark_1D_params()
  @timeit to "Warm-up 1D test" begin
    @test DamBreakSW.main(params_1D) == nothing
  end
  params_1D = DamBreakSW.DamBreak_benchmark_1D_params(
    nₓ=400,
    verbose=true
  )
  @timeit to "Fine 1D test" begin
    @test DamBreakSW.main(params_1D) == nothing
  end
end

show(to; compact=true)
