using DrWatson
@quickactivate "DamBreakSW"

include(srcdir("DamBreakSW.jl"))

params = DamBreakSW.Mesh_params(
  H=0.15,
  h=0.015,
  decay_factor_left=3,
  decay_factor_right=2.5,
  dxLeft=3.0,
  dyTop=0.8,
  dyBottom=0.8
)
DamBreakSW.generate_mesh(params)
