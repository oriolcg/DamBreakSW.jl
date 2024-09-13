"""
    struct physics_params

A structure to store the parameters for physical quantities.
The parameters are:
  - `ν::Float64`: Viscosity
  - `Cd::Float64`: Drag coefficient
  - `g::Float64`: Gravity
"""
@with_kw struct physics_params
  h₀⬆::Float64 = 10.0   # Initial height of the water upstream
  h₀⬇::Float64 = 5.0    # Initial height of the water downstream
  ν::Float64 = 1.0      # Viscosity
  Cd::Float64 = 0.0     # Drag coefficient
  g::Float64 = 9.8      # Gravity
end

# FE Spaces
function get_FESpaces(Ω,D,order::Int,DTags,DMasks,DValues,::Union{Val{:Galerkin},Val{:ASGS},Val{:Smagorinsky}})
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  refFEₕ = ReferenceFE(lagrangian,Float64,order-1)
  Vᵤ = TestFESpace(Ω,refFEᵤ,dirichlet_tags=DTags,dirichlet_masks=DMasks)
  Vₕ = TestFESpace(Ω,refFEₕ;conformity=:H1)
  Uᵤ = TransientTrialFESpace(Vᵤ,DValues)
  Uₕ = TransientTrialFESpace(Vₕ)
  Y = MultiFieldFESpace([Vᵤ,Vₕ])
  X = TransientMultiFieldFESpace([Uᵤ,Uₕ])
  return X,Y
end

function get_FESpaces(Ω,D,order::Int,DTags,DMasks,DValues,::Val{:VectorInvariant})
  refFEᵤ = ReferenceFE(raviart_thomas,Float64,order-1)
  refFEₕ = ReferenceFE(lagrangian,Float64,order-1)
  refFEq = ReferenceFE(lagrangian,Float64,order)
  Vᵤ = TestFESpace(Ω,refFEᵤ,conformity=:Hdiv)#,dirichlet_tags=DTags,dirichlet_masks=DMasks)
  Vₕ = TestFESpace(Ω,refFEₕ;conformity=:L2)
  Vq = TestFESpace(Ω,refFEq;conformity=:H1)
  Uᵤ = TransientTrialFESpace(Vᵤ)#,DValues)
  Uₕ = TransientTrialFESpace(Vₕ)
  Uq = TransientTrialFESpace(Vq)
  Y = MultiFieldFESpace([Vᵤ,Vₕ,Vᵤ,Vₕ,Vq])
  X = TransientMultiFieldFESpace([Uᵤ,Uₕ,Uᵤ,Uₕ,Uq])
  return (X,Y)
end

# Residual functions
"""
    get_forms(::Val{:Galerkin},params::physics_params)

Get the FE operator forms for the standard Galerkin method.
"""
function get_forms(measures,normals,D,::Val{:Galerkin},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  absᵤ(u) = √(u⋅u + 1.0e-8)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  strs(∇u,∇v) = ( ν*(∇u+∇u') - 2/3*ν*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(v⋅∇h)
  # grad(h,∇v) = -g*h*(tr(∇v))
  # gradΓ(h,v,n) = g*h*(v⋅n)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w
  # convₕ(u,h,∇u,∇h,w) = ((h+h₀⬇)*u)⋅∇w

  # Residual form
  dΩ,dΓwall, = measures
  nwall, = normals
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        (strs∘(∇(u),∇(v))) +
                        (drag∘(u,h,v)) +
                        # (grad∘(h,∇(v))) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ# +
                    #  ∫( (gradΓ∘(h,v,nwall)) )dΓwall
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))

  return m,a,res

end

"""
    get_forms(::Val{:ASGS},params::physics_params)

Get the FE operator forms for the stabilized formulation using Algebraic Subgrid Scales (ASGS).
"""
function get_forms(measures,normals,D,::Val{:ASGS},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  absᵤ(u) = (u⋅u + 1.0e-8).^(1/2)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  strs(∇u,∇v) = ( ν*(∇u+∇u') - 2/3*ν*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(∇h⋅v)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w

  # Stabilization
  c₁ = 12.0; c₂ = 2.0; c₃ = 1.0
  Rₕ(u,h,hₜ,∇u,∇h) = hₜ + u⋅∇h + (h+h₀⬇)*(tr(∇u))
  Rᵤ(u,h,uₜ,∇u,∇h) = uₜ + u⋅∇u + g*∇h + Cd/(h+h₀⬇)*(absᵤ(u))*u
  Lᵤᵃ(u,∇v,∇w) = - u⋅∇v - g*∇w
  Lₕᵃ(u,h,∇v,∇w) = - u⋅∇w - (h+h₀⬇)*(tr(∇v))
  τᵤ(a,h,Δxₒ) = 1.0 / ((c₁*ν / (Δxₒ^2)) + (c₂*absᵤ(a) / Δxₒ) + (c₃*Cd*g*absᵤ(a) / (h+1.0e-8)))
  τₕ(a,h,Δxₒ) = (Δxₒ^2)/(c₁*τᵤ(a,h,Δxₒ))
  stabₕ(u,h,hₜ,∇u,∇h,∇v,∇w,Δxₒ) = (τₕ∘(u,h,Δxₒ))*((Rₕ∘(u,h,hₜ,∇u,∇h))*Lₕᵃ(u,h,∇v,∇w))
  stabᵤ(u,h,uₜ,∇u,∇h,∇v,∇w,Δxₒ) = (τᵤ∘(u,h,Δxₒ))*((Rᵤ∘(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w))

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = get_triangulation(dΩ.quad)
  Δxₒ = lazy_map(dx->dx^(1/D),get_cell_measure(Ω))

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        (strs∘(∇(u),∇(v))) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) -
                        (stabᵤ(u,h,∂t(u),∇(u),∇(h),∇(v),∇(w),Δxₒ)) )dΩ#-
                        # (stabₕ(u,h,∂t(h),∇(u),∇(h),∇(v),∇(w),Δxₒ)) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))

  return m,a,res

end

"""
    get_forms(::Val{:Smagorinsky},params::physics_params)

Get the operator forms for the Smagorinsky method.
"""
function get_forms(measures,normals,D,::Val{:Smagorinsky},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  cₛ = 0.164
  νₜ(εᵤ,Δx₀) = (cₛ*Δx₀)^2*(√(2*(εᵤ⊙εᵤ)+1.0e-8 ))
  absᵤ(u) = √(u⋅u + 1.0e-8)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  strs(∇u,∇v,εᵤ,Δx₀) = ( (ν+νₜ(εᵤ,Δx₀))*(∇u+∇u') - 2/3*(ν+νₜ(εᵤ,Δx₀))*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(v⋅∇h)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = get_triangulation(dΩ.quad)
  Δx₀ = lazy_map(dx->dx^(1/D),get_cell_measure(Ω))

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        (strs∘(∇(u),∇(v),ε(u),Δx₀)) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))

  return m,a,res

end

"""
    get_forms(::Val{:Smagorinsky},params::physics_params)

Get the operator forms for the Smagorinsky method.
"""
function get_forms(measures::Tuple{Vararg{GridapDistributed.DistributedMeasure}},normals,D,::Val{:Smagorinsky},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  cₛ = 0.164
  νₜ(εᵤ,Δx₀) = (cₛ*Δx₀)*(√(2*(εᵤ⊙εᵤ)+1.0e-8 ))
  absᵤ(u) = (u⋅u + 1.0e-8).^(1/2)
  dabsᵤ(u,du) = 1/absᵤ(u)*(du⋅u)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  strs(∇u,∇v,εᵤ,Δx₀) = ( (ν+νₜ(εᵤ,Δx₀))*(∇u+∇u') - 2/3*(ν+νₜ(εᵤ,Δx₀))*tr(∇u)*I) ⊙ ∇v
  strs(∇u,∇v,εᵤ) = ( (ν+νₜ(εᵤ,0.01))*(∇u+∇u') - 2/3*(ν+νₜ(εᵤ,0.01))*tr(∇u)*I) ⊙ ∇v
  strs(∇u,∇v) = ( (ν)*(∇u+∇u') - 2/3*(ν)*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  dudrag(u,du,h,v) = Cd/(h+h₀⬇)*((dabsᵤ(u,du))*(u⋅v)+(absᵤ(u))*(du⋅v))
  dhdrag(u,h,dh,v) = -Cd/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*(u⋅v)*dh
  grad(∇h,v) = g*(v⋅∇h)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w
  dhconvₕ(u,dh,∇u,∇dh,w) = ((u⋅∇dh) + (dh)*tr(∇u))*w

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = dΩ.trian
  Δx₀ = map(Ω.trians) do trian
    lazy_map(dx->dx^(1/D),get_cell_measure(trian))
  end

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        # (strs∘(∇(u),∇(v),ε(u),Δx₀)) +
                        (strs∘(∇(u),∇(v),ε(u))) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))
  jac(t,(u,h),(du,dh),(v,w)) =
    ∫( (convᵤ(du,∇(u),v))  )dΩ +
    ∫( (convᵤ(u,∇(du),v)) )dΩ +
    ∫( (strs(∇(du),∇(v))) )dΩ +
    ∫( (dudrag(u,du,h,v)) )dΩ +
    ∫( (dhdrag(u,h,dh,v)) )dΩ +
    ∫( (grad(∇(dh),v)) )dΩ +
    ∫( (convₕ(du,h,∇(du),∇(h),w)) )dΩ +
    ∫( (dhconvₕ(u,dh,∇(u),∇(dh),w)) )dΩ
  jac_t(t,(uₜ,hₜ),(duₜ,dhₜ),(v,w)) = m(t,(duₜ,dhₜ),(v,w))

  return m,a,res,jac,jac_t

end

"""
    get_forms(::Val{:VectorInvariant},params::physics_params)

Get the operator forms for the vector-invariant formulation.
"""
function get_forms(measures,normals,D,::Val{:VectorInvariant},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  vecPerp(u) = VectorValue(-deepcopy(u[2]),deepcopy(u[1]))
  gradPerp(∇ϕ::VectorValue{2}) = VectorValue( -deepcopy(∇ϕ[2]), deepcopy(∇ϕ[1]))

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = get_triangulation(dΩ.quad)
  Δx₀ = lazy_map(dx->dx^(1/D),get_cell_measure(Ω))

  # Residual form
  @unpack Δt = ode_solver_params
  m(t,(uₜ,hₜ,),(v,w,)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h,F,Φ,q),(v,w,s,ψ,p)) = (
    ∫( (∇⋅F)*w  )dΩ
  - ∫( (∇⋅v)*Φ  )dΩ
  + ∫( (q - 0.5*Δt*(u⋅∇(q)) )*(vecPerp∘(F)⋅v)  )dΩ
  + ∫( F⋅s - h*(u⋅s) )dΩ
  + ∫( Φ*ψ - (0.5*(u⋅u) + g*h)*ψ )dΩ
  + ∫( q*h*p + gradPerp∘(∇(p))⋅u )dΩ
  )
  res(t,(u,h,F,Φ,q),(v,w,s,ψ,p)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))

  return m,a,res

end
