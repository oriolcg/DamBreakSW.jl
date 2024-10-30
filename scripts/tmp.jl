module TMP

using Gridap
using GridapGmsh
using Gridap.TensorValues: meas
using Gridap.FESpaces: interpolate_everywhere!

function main(file,nT)

  𝒯 = GmshDiscreteModel(file)

  # Triangulations
  Ω = Interior(𝒯)
  Γₒᵤₜ = Boundary(𝒯,tags="outflow")
  Γw = Boundary(𝒯,tags="walls")
  Γₛ = Boundary(𝒯,tags="sides")
  Λ = Skeleton(Ω)

  # Parameters & Boundaries
  g = 9.81
  H=0.03
  h₀(x,t) = (0.63-H) * (x[1]<=0)
  h₀(t::Real) = x->h₀(x,t)
  u₀(x,t) = VectorValue(0.0,0.0)
  u₀(t::Real) = x->u₀(x,t)

  ν = 1.0e-6 #3.225
  cD = 0.0127
  g = 9.81
  hₒᵤₜ(x,t) = 0.001
  hₒᵤₜ(t) = x -> hₒᵤₜ(x,t)
  I = TensorValue(1.0,0.0,0.0,1.0)

  # FE Spaces
  order = 1
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  refFEₕ = ReferenceFE(lagrangian,Float64,order-1)
  Vᵤ = TestFESpace(Ω,refFEᵤ,dirichlet_tags="walls")
  Vₕ = TestFESpace(Ω,refFEₕ;conformity=:L2)
  Uᵤ = TransientTrialFESpace(Vᵤ,u₀)
  Uₕ = TransientTrialFESpace(Vₕ)
  Y = MultiFieldFESpace([Vᵤ,Vₕ])
  X = TransientMultiFieldFESpace([Uᵤ,Uₕ])

  # Measure (quadrature rules for numeric integration)
  dΩ = Measure(Ω,2*order)
  dΓₒᵤₜ = Measure(Γₒᵤₜ,2*order)
  dΛ = Measure(Λ,2*order)
  dΓw = Measure(Γw,2*order)
  dΓₛ = Measure(Γₛ,2*order)

  # Normal vectors to the boundaries
  nΓₒᵤₜ = get_normal_vector(Γₒᵤₜ)
  nΛ = get_normal_vector(Λ)
  nΓw = get_normal_vector(Γw)
  nΓₛ = get_normal_vector(Γₛ)

  # Mesh size
  Δxₒ = CellField(lazy_map(dx->dx^(1/2),get_cell_measure(Ω)),Ω)
  ΔxΛ = 0.2#get_cell_measure(Λ)

  #Max function of water elevation
  maxh(h)=max(h+H,1.0e-10)

  # Stabilization
  c₁ = 12.0; c₂ = 2.0; c₃ = 1.0
  global uₙ = interpolate_everywhere(u₀(0.0),Uᵤ(0.0))
  global hₙ = interpolate_everywhere(h₀(0.0),Uₕ(0.0))
  #meas∘(u) = √(u⋅u)#meas∘(u)
  dmeasu(u,du) = (u⋅du)/(√(u⋅u)+1e-14)
  Rₕ(u,h) = ∂t(h) + ∇(h)'⋅u + (h+H)*(∇⋅u)
  Rᵤ(u,h) = ∂t(u) + ∇(u)'⋅u + g*∇(h) + cD/(h+H)*u*(meas∘(u))
  dRₕ(u,h,du,dh) = ∇(dh)⋅u + dh*(∇⋅u) + ∇(h)⋅du + (h+H)*(∇⋅du)
  dRᵤ(u,h,du,dh) = ∇(du)'⋅u + + ∇(u)'⋅du + g*∇(dh) + cD/(h+H)*du*(meas∘(u)) + cD/(h+H)*u*(dmeasu∘(u,du)) - cD/(h*h+1e-14)*u*dh*(meas∘(u))
  Lᵤᵃ(v,w) = - ∇(v)'⋅uₙ - g*∇(w)
  Lₕᵃ(v,w) = - ∇(w)⋅uₙ - H*(∇⋅v)
  τᵤ(a,h) = 1.0 / (c₁*ν/(Δxₒ*Δxₒ) + c₂*a/Δxₒ + c₃*cD*g*a/(h+1e-14))
  τₕ(a,h) = (Δxₒ*Δxₒ)/(c₁*τᵤ(a,h))
  dτᵤdu(a,h,da) = - τᵤ(a,h)*τᵤ(a,h) * (c₂/Δxₒ + c₃*cD*g/(h+1e-14))*da
  dτᵤdh(a,h,dh) = τᵤ(a,h)*τᵤ(a,h) * c₃*cD*g*a/(h*h+1e-14)*dh
  dτₕdu(a,h,da) = τₕ(a,h)/τᵤ(a,h)*dτᵤdu(a,h,da)
  dτₕdh(a,h,dh) = τₕ(a,h)/τᵤ(a,h)*dτᵤdh(a,h,dh)
  γ = 1.0/ΔxΛ

  # Weak form
  # =========
  res(t,(u,h),(v,w)) = ∫( ∂t(h)*w - (h+H)*u⋅∇(w) +
      (∂t(u) + ∇(u)'⋅u + cD*((meas∘(u))/(h+H))*u) ⋅ v - g*h*(∇⋅v) +
      ν*( (∇(u)+∇(u)') - 2/3*(∇⋅u)*I ) ⊙ ∇(v) -
      Rₕ(u,h) * ((τₕ(meas∘u,h))*Lₕᵃ(v,w)) -
      Rᵤ(u,h) ⋅ ((τᵤ(meas∘u,h))*Lᵤᵃ(v,w)) )dΩ +
    ∫( g*(hₒᵤₜ(t)*(v⋅nΓₒᵤₜ)) + hₒᵤₜ(t)*(u⋅nΓₒᵤₜ)*w + H*(u⋅nΓₒᵤₜ)*w )dΓₒᵤₜ +
    #  ∫( g*h*(v⋅nΓw) )dΓw +
    ∫( g*h*(v⋅nΓₛ) )dΓₛ+
    ∫( mean((h+H)*u)⋅jump(w*nΛ) + γ*jump(h*nΛ)⋅jump(w*nΛ) +
      mean(g*h)*jump(v⋅nΛ) )dΛ
  jac(t,(u,h),(du,dh),(v,w)) = ∫(
      - ((h+H)*du + dh*u)⋅∇(w) +
      (∇(du)'⋅u + ∇(u)'⋅du + cD*((meas∘(u))/(h+H))*du + cD/(h+H)*u*(dmeasu∘(u,du)) - cD*((meas∘(u))/(h+H))*u*dh ) ⋅ v - g*dh*(∇⋅v) +
      ν*( (∇(du)+∇(du)') - 2/3*(∇⋅du)*I ) ⊙ ∇(v) -
      dRₕ(u,h,du,dh) * ((τₕ(meas∘u,h))*Lₕᵃ(v,w)) -
      dRᵤ(u,h,du,dh) ⋅ ((τᵤ(meas∘u,h))*Lᵤᵃ(v,w)) -
      Rₕ(u,h) * ((dτₕdu(meas∘u,h,dmeasu∘(u,du)) + dτₕdh(meas∘u,h,dh))*Lₕᵃ(v,w)) -
      Rᵤ(u,h) ⋅ ((dτᵤdu(meas∘u,h,dmeasu∘(u,du)) + dτᵤdh(meas∘u,h,dh))*Lᵤᵃ(v,w)) )dΩ +
    ∫( hₒᵤₜ(t)*(du⋅nΓₒᵤₜ)*w + H*(du⋅nΓₒᵤₜ)*w )dΓₒᵤₜ +
    #  ∫( g*dh*(v⋅nΓw) )dΓw +
    ∫( g*dh*(v⋅nΓₛ) )dΓₛ +
    ∫( mean(dh*u)⋅jump(w*nΛ) + mean((h+H)*du)⋅jump(w*nΛ) +
      γ*jump(dh*nΛ)⋅jump(w*nΛ) +
      mean(g*dh)*jump(v⋅nΛ) )dΛ
  jac_t(t,(u,h),(dut,dht),(v,w)) = ∫(
    dht * w +
    dut ⋅ v -
    dht * (τₕ(meas∘u,h)*Lₕᵃ(v,w)) -
    dut ⋅ (τᵤ(meas∘u,h)*Lᵤᵃ(v,w)) )dΩ

    op = TransientFEOperator(res,(jac,jac_t),X,Y)

  # Solver
  nls = NLSolver(show_trace=true,iterations=10,ftol=1.0e-6)
  t₀ = 0.0
  Δt= 0.05
  T = nT*Δt
  # ode_scheme = ThetaMethod(nls, Δt, 0.5)
  ode_scheme = GeneralizedAlpha1(nls,Δt,0.0)
  xₕ₀ = interpolate_everywhere([u₀(0.0),h₀(0.0)],X(0.0))
  xdotₕ₀ = interpolate_everywhere([VectorValue(0.0,0.0),0.0],X(0.0))
  xₕₜ = solve(ode_scheme,op,t₀,T,(xₕ₀,xdotₕ₀))

  output_files = createpvd("data\\sims\\tmp\\tmp", append=false) do pvd
    for (t,xₕ) in xₕₜ
      println("Time: $t")
      uₕ, hₕ = xₕ
      pvd[t] = createvtk(Ω,"tmp_$t.vtu",cellfields=["u"=>uₕ,"h"=>hₕ])
      global uₙ, hₙ
      interpolate_everywhere!(uₙ,get_free_dof_values(uₕ),get_dirichlet_dof_values(Uᵤ(t)),Uᵤ(t))
      interpolate_everywhere!(hₙ,get_free_dof_values(hₕ),get_dirichlet_dof_values(Uₕ(t)),Uₕ(t))
    end
  end
end

main("data\\meshes\\DamBreak_building_0.0_coarse.msh",1)
main("data\\meshes\\DamBreak_building_0.0_medium.msh",20)

end
