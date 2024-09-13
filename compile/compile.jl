using PackageCompiler
create_sysimage(:DamBreakSW,
  sysimage_path=joinpath(@__DIR__,"..","DamBreakSW.so"),
  precompile_execution_file=joinpath(@__DIR__,"warmup.jl"))
