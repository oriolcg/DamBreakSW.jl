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
function get_FESpaces(Ω,D,order::Int,DTags,DMasks,DValues,::Union{Val{:Galerkin},Val{:Smagorinsky}})
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  refFEₕ = ReferenceFE(lagrangian,Float64,order-1)
  Vₕ = TestFESpace(Ω,refFEₕ;conformity=:H1)
  if isempty(DValues)
    Vᵤ = TestFESpace(Ω,refFEᵤ)
    Uᵤ = TransientTrialFESpace(Vᵤ)
  else
    Vᵤ = TestFESpace(Ω,refFEᵤ,dirichlet_tags=DTags,dirichlet_masks=DMasks)
    Uᵤ = TransientTrialFESpace(Vᵤ,DValues)
  end
  Uₕ = TransientTrialFESpace(Vₕ)
  Y = MultiFieldFESpace([Vᵤ,Vₕ])
  X = TransientMultiFieldFESpace([Uᵤ,Uₕ])
  return X,Y
end

function get_FESpaces(Ω,D,order::Int,DTags,DMasks,DValues,::Val{:ASGS})
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  refFEₕ = ReferenceFE(lagrangian,Float64,order)
  Vₕ = TestFESpace(Ω,refFEₕ;conformity=:H1)
  if isempty(DValues)
    Vᵤ = TestFESpace(Ω,refFEᵤ)
    Uᵤ = TransientTrialFESpace(Vᵤ)
  else
    Vᵤ = TestFESpace(Ω,refFEᵤ,dirichlet_tags=DTags,dirichlet_masks=DMasks)
    Uᵤ = TransientTrialFESpace(Vᵤ,DValues)
  end
  Uₕ = TransientTrialFESpace(Vₕ)
  Y = MultiFieldFESpace([Vᵤ,Vₕ])
  X = TransientMultiFieldFESpace([Uᵤ,Uₕ])
  return X,Y
end

function get_FESpaces(Ω,D,order::Int,DTags,DMasks,DValues,::Val{:ASGS_L2})
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  refFEₕ = ReferenceFE(lagrangian,Float64,order-1)
  Vₕ = TestFESpace(Ω,refFEₕ;conformity=:L2)
  if isempty(DValues)
    Vᵤ = TestFESpace(Ω,refFEᵤ)
    Uᵤ = TransientTrialFESpace(Vᵤ)
  else
    Vᵤ = TestFESpace(Ω,refFEᵤ,dirichlet_tags=DTags,dirichlet_masks=DMasks)
    Uᵤ = TransientTrialFESpace(Vᵤ,DValues)
  end
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

function get_FESpaces(Ω,D,order::Int,DTags,DMasks,DValues,::Val{:conservative_Galerkin})
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{D+1,Float64},order)
  if isempty(DValues)
    V = TestFESpace(Ω,refFEᵤ)
    U = TransientTrialFESpace(V)
  else
    V = TestFESpace(Ω,refFEᵤ,dirichlet_tags=DTags,dirichlet_masks=DMasks)
    U = TransientTrialFESpace(V,DValues)
  end
  return U,V
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
  @unpack Δt = ode_solver_params
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
  dudrag(u,du,h,v) = Cd/(h+h₀⬇)*((duabsᵤ(u,du))*(u⋅v)+(absᵤ(u))*(du⋅v))
  dhdrag(u,h,dh,v) = -Cd/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*(u⋅v)*dh
  dhconvₕ(u,dh,∇u,∇dh,w) = ((u⋅∇dh) + (dh)*tr(∇u))*w
  duabsᵤ(u,du) = 1.0 / (u⋅u + 1.0e-8).^(1/2) * (u⋅du)

  # Stabilization
  c₁ = 12.0; c₂ = 2.0; c₃ = 1.0
  Rₕ(u,h,hₜ,∇u,∇h) = hₜ + u⋅∇h + (h+h₀⬇)*(tr(∇u))
  Rᵤ(u,h,uₜ,∇u,∇h) = uₜ + u⋅∇u + g*∇h + Cd/(h+h₀⬇)*(absᵤ(u))*u
  Lᵤᵃ(u,∇v,∇w) = - u⋅∇v - g*∇w
  Lₕᵃ(u,h,∇v,∇w) = - u⋅∇w - (h+h₀⬇)*(tr(∇v))
  τᵤinv(a,h,Δx₀) = 2/Δt + (c₁*ν / (Δx₀*Δx₀)) + (c₂*absᵤ(a) / Δx₀) + (c₃*Cd*g*absᵤ(a) / (h+1.0e-8))
  τᵤ(a,h,Δx₀) = 1.0 / τᵤinv(a,h,Δx₀)
  # τᵤ₂ = Δt/2
  # τᵤ(a,h,Δx₀) = 1.0 / (1.0/τᵤinv(a,h,Δx₀) + 1.0/τᵤ₂)
  τₕ(a,h,Δx₀) = τᵤ(a,h,Δx₀)#(Δx₀^2)/(c₁*τᵤ(a,h,Δx₀))
  stabₕ(u,h,hₜ,∇u,∇h,∇v,∇w,Δx₀) = (τₕ∘(u,h,Δx₀))*((Rₕ∘(u,h,hₜ,∇u,∇h))*Lₕᵃ(u,h,∇v,∇w))
  stabᵤ(u,h,uₜ,∇u,∇h,∇v,∇w,Δx₀) = (τᵤ∘(u,h,Δx₀))*((Rᵤ∘(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w))
  # dustabᵤ(u,du,h,uₜ,∇u,∇du,∇h,∇v,∇w,Δx₀) =
  #   (duτᵤ(u,du,h,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
  #   (τᵤ(u,h,Δx₀))*((duRᵤ(u,du,h,uₜ,∇u,∇du,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
  #   (τᵤ(u,h,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅duLᵤᵃ(du,∇v,∇w))
  # dhstabᵤ(u,h,dh,uₜ,∇u,∇h,∇dh,∇v,∇w,Δx₀) =
  #   (dhτᵤ(u,h,dh,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
  #   (τᵤ(u,h,Δx₀))*((dhRᵤ(u,h,dh,uₜ,∇u,∇dh))⋅Lᵤᵃ(u,∇v,∇w))
  # duₜstabᵤ(u,h,duₜ,∇u,∇h,∇v,∇w,Δx₀) = (τᵤ(u,h,Δx₀))*((duₜRᵤ(duₜ))⋅Lᵤᵃ(u,∇v,∇w))
  # duτᵤ(a,da,h,Δx₀) = -1.0 / (τᵤinv(a,h,Δx₀)*τᵤinv(a,h,Δx₀)) * duτᵤinv(a,da,h,Δx₀)
  # duτᵤinv(a,da,h,Δx₀) = (c₂*duabsᵤ(a,da) / Δx₀) + (c₃*Cd*g*duabsᵤ(a,da) / (h+1.0e-8))
  # dhτᵤ(a,h,dh,Δx₀) = -1.0 / (τᵤinv(a,h,Δx₀)*τᵤinv(a,h,Δx₀)) * dhτᵤinv(a,h,dh,Δx₀)
  # dhτᵤinv(a,h,dh,Δx₀) = -1.0*(c₃*Cd*g*absᵤ(a) / ((h+1.0e-8)*(h+1.0e-8))) * dh
  # duRᵤ(u,du,h,uₜ,∇u,∇du,∇h) = u⋅∇du + du⋅∇u + Cd/(h+h₀⬇)*(duabsᵤ(u,du))*u + Cd/(h+h₀⬇)*(absᵤ(u))*du
  # duLᵤᵃ(du,∇v,∇w) = - du⋅∇v
  # dhRᵤ(u,h,dh,uₜ,∇u,∇dh) = g*∇dh - Cd/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*u *dh
  # duₜRᵤ(duₜ) = duₜ

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = get_triangulation(dΩ.quad)
  Δx₀ = CellField(lazy_map(dx->dx^(1/D),get_cell_measure(Ω)),Ω)

  # Residual form
  res(t,(u,h),(v,w)) = ∫( ∂t(u)⋅v + ∂t(h)*w )dΩ +
                       ∫( (convᵤ∘(u,∇(u),v)) )dΩ +
                       ∫( (strs∘(∇(u),∇(v))) )dΩ +
                       ∫( (drag∘(u,h,v)) )dΩ +
                       ∫( (grad∘(∇(h),v)) )dΩ +
                       ∫( (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ -
                       ∫( (stabᵤ(u,h,∂t(u),∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ -
                       ∫( (stabₕ(u,h,∂t(h),∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ
  # jac(t,(u,h),(du,dh),(v,w)) =
  #                      ∫( (convᵤ(du,∇(u),v))  )dΩ +
  #                      ∫( (convᵤ(u,∇(du),v)) )dΩ +
  #                      ∫( (strs(∇(du),∇(v))) )dΩ +
  #                      ∫( (dudrag(u,du,h,v)) )dΩ +
  #                      ∫( (dhdrag(u,h,dh,v)) )dΩ +
  #                      ∫( (grad(∇(dh),v)) )dΩ +
  #                      ∫( (convₕ(du,h,∇(du),∇(h),w)) )dΩ +
  #                      ∫( (dhconvₕ(u,dh,∇(u),∇(dh),w)) )dΩ -
  #                      ∫( (dustabᵤ(u,du,h,∂t(u),∇(u),∇(du),∇(h),∇(v),∇(w),Δx₀)) )dΩ -
  #                      ∫( (dhstabᵤ(u,h,dh,∂t(u),∇(u),∇(h),∇(dh),∇(v),∇(w),Δx₀)) )dΩ
  # jac_t(t,(u,h),(duₜ,dhₜ),(v,w)) =
  #   ∫( duₜ⋅v + dhₜ*w )dΩ -
  #   ∫( (duₜstabᵤ(u,h,duₜ,∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ

  return nothing,nothing,res#,(jac,jac_t)

end

"""
    get_forms(::Val{:ASGS_L2},params::physics_params)

Get the FE operator forms for the stabilized formulation using Algebraic Subgrid Scales (ASGS) with L2-conforming spaces in h.
"""
function get_forms(measures,normals,D,::Val{:ASGS},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  @unpack Δt = ode_solver_params
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
  dudrag(u,du,h,v) = Cd/(h+h₀⬇)*((duabsᵤ(u,du))*(u⋅v)+(absᵤ(u))*(du⋅v))
  dhdrag(u,h,dh,v) = -Cd/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*(u⋅v)*dh
  dhconvₕ(u,dh,∇u,∇dh,w) = ((u⋅∇dh) + (dh)*tr(∇u))*w
  duabsᵤ(u,du) = 1.0 / (u⋅u + 1.0e-8).^(1/2) * (u⋅du)

  # Stabilization
  c₁ = 12.0; c₂ = 2.0; c₃ = 1.0
  Rₕ(u,h,hₜ,∇u,∇h) = hₜ + u⋅∇h + (h+h₀⬇)*(tr(∇u))
  Rᵤ(u,h,uₜ,∇u,∇h) = uₜ + u⋅∇u + g*∇h + Cd/(h+h₀⬇)*(absᵤ(u))*u
  Lᵤᵃ(u,∇v,∇w) = - u⋅∇v - g*∇w
  Lₕᵃ(u,h,∇v,∇w) = - u⋅∇w - (h+h₀⬇)*(tr(∇v))
  τᵤinv(a,h,Δx₀) = 2/Δt + (c₁*ν / (Δx₀*Δx₀)) + (c₂*absᵤ(a) / Δx₀) + (c₃*Cd*g*absᵤ(a) / (h+1.0e-8))
  τᵤ(a,h,Δx₀) = 1.0 / τᵤinv(a,h,Δx₀)
  # τᵤ₂ = Δt/2
  # τᵤ(a,h,Δx₀) = 1.0 / (1.0/τᵤinv(a,h,Δx₀) + 1.0/τᵤ₂)
  τₕ(a,h,Δx₀) = τᵤ(a,h,Δx₀)#(Δx₀^2)/(c₁*τᵤ(a,h,Δx₀))
  stabₕ(u,h,hₜ,∇u,∇h,∇v,∇w,Δx₀) = (τₕ∘(u,h,Δx₀))*((Rₕ∘(u,h,hₜ,∇u,∇h))*Lₕᵃ(u,h,∇v,∇w))
  stabᵤ(u,h,uₜ,∇u,∇h,∇v,∇w,Δx₀) = (τᵤ∘(u,h,Δx₀))*((Rᵤ∘(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w))
  # dustabᵤ(u,du,h,uₜ,∇u,∇du,∇h,∇v,∇w,Δx₀) =
  #   (duτᵤ(u,du,h,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
  #   (τᵤ(u,h,Δx₀))*((duRᵤ(u,du,h,uₜ,∇u,∇du,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
  #   (τᵤ(u,h,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅duLᵤᵃ(du,∇v,∇w))
  # dhstabᵤ(u,h,dh,uₜ,∇u,∇h,∇dh,∇v,∇w,Δx₀) =
  #   (dhτᵤ(u,h,dh,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
  #   (τᵤ(u,h,Δx₀))*((dhRᵤ(u,h,dh,uₜ,∇u,∇dh))⋅Lᵤᵃ(u,∇v,∇w))
  # duₜstabᵤ(u,h,duₜ,∇u,∇h,∇v,∇w,Δx₀) = (τᵤ(u,h,Δx₀))*((duₜRᵤ(duₜ))⋅Lᵤᵃ(u,∇v,∇w))
  # duτᵤ(a,da,h,Δx₀) = -1.0 / (τᵤinv(a,h,Δx₀)*τᵤinv(a,h,Δx₀)) * duτᵤinv(a,da,h,Δx₀)
  # duτᵤinv(a,da,h,Δx₀) = (c₂*duabsᵤ(a,da) / Δx₀) + (c₃*Cd*g*duabsᵤ(a,da) / (h+1.0e-8))
  # dhτᵤ(a,h,dh,Δx₀) = -1.0 / (τᵤinv(a,h,Δx₀)*τᵤinv(a,h,Δx₀)) * dhτᵤinv(a,h,dh,Δx₀)
  # dhτᵤinv(a,h,dh,Δx₀) = -1.0*(c₃*Cd*g*absᵤ(a) / ((h+1.0e-8)*(h+1.0e-8))) * dh
  # duRᵤ(u,du,h,uₜ,∇u,∇du,∇h) = u⋅∇du + du⋅∇u + Cd/(h+h₀⬇)*(duabsᵤ(u,du))*u + Cd/(h+h₀⬇)*(absᵤ(u))*du
  # duLᵤᵃ(du,∇v,∇w) = - du⋅∇v
  # dhRᵤ(u,h,dh,uₜ,∇u,∇dh) = g*∇dh - Cd/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*u *dh
  # duₜRᵤ(duₜ) = duₜ

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = get_triangulation(dΩ.quad)
  Δx₀ = CellField(lazy_map(dx->dx^(1/D),get_cell_measure(Ω)),Ω)

  # ∂ₜh + ∇⋅(u*(h+h₀⬇)) = 0
  # ∂ₜu + u⋅∇u - ∇⋅(2με(u)) + g∇h + + Cd/(h+h₀⬇)*(absᵤ(u))*u = 0

  # (∂ₜh,w) + (∇⋅(u*(h+h₀⬇)),w) = 0
  # (∂ₜu,v) + (u⋅∇u,v) - (∇⋅(2με(u)),v) + (g∇h,v) +  (Cd/(h+h₀⬇)*(absᵤ(u))*u,v) = 0

  # (∂ₜh,w) - ((u*(h+h₀⬇)), ∇w) + ((u*(h+h₀⬇))⋅n,w)_Γ = 0
  # (∂ₜu,v) + (u⋅∇u,v) + (2με(u)),∇v) - (2με(u)⋅n,v)_Γ - (gh,∇⋅v) + (gh*n,v)_Γ +  (Cd/(h+h₀⬇)*(absᵤ(u))*u,v) = 0

  # Residual form
  res(t,(u,h),(v,w)) = ∫( ∂t(u)⋅v + ∂t(h)*w )dΩ +
                       ∫( (convᵤ∘(u,∇(u),v)) )dΩ +
                       ∫( (strs∘(∇(u),∇(v))) )dΩ +
                       ∫( (drag∘(u,h,v)) )dΩ -
                       ∫( g*h*(∇⋅v) )dΩ +
                       # BC for h go here (v=0 on Dirichlet u)-
                       ∫( (u*(h+h₀⬇))⋅(∇(w)) )dΩ +
                       ∫( jump((u*(h+h₀⬇))⋅nΛ)*mean(w) )dΛ +
                       # BC for convective term goes here (w=0 on Dirichlet h) -
                      #  ∫( (stabᵤ(u,h,∂t(u),∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ -
                      #  ∫( (stabₕ(u,h,∂t(h),∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ
  # jac(t,(u,h),(du,dh),(v,w)) =
  #                      ∫( (convᵤ(du,∇(u),v))  )dΩ +
  #                      ∫( (convᵤ(u,∇(du),v)) )dΩ +
  #                      ∫( (strs(∇(du),∇(v))) )dΩ +
  #                      ∫( (dudrag(u,du,h,v)) )dΩ +
  #                      ∫( (dhdrag(u,h,dh,v)) )dΩ +
  #                      ∫( (grad(∇(dh),v)) )dΩ +
  #                      ∫( (convₕ(du,h,∇(du),∇(h),w)) )dΩ +
  #                      ∫( (dhconvₕ(u,dh,∇(u),∇(dh),w)) )dΩ -
  #                      ∫( (dustabᵤ(u,du,h,∂t(u),∇(u),∇(du),∇(h),∇(v),∇(w),Δx₀)) )dΩ -
  #                      ∫( (dhstabᵤ(u,h,dh,∂t(u),∇(u),∇(h),∇(dh),∇(v),∇(w),Δx₀)) )dΩ
  # jac_t(t,(u,h),(duₜ,dhₜ),(v,w)) =
  #   ∫( duₜ⋅v + dhₜ*w )dΩ -
  #   ∫( (duₜstabᵤ(u,h,duₜ,∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ

  return nothing,nothing,res#,(jac,jac_t)

end

function get_forms(measures::Tuple{Vararg{GridapDistributed.DistributedMeasure}},normals,D,::Val{:ASGS},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  @unpack Δt = ode_solver_params
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
  dudrag(u,du,h,v) = Cd/(h+h₀⬇)*((duabsᵤ(u,du))*(u⋅v)+(absᵤ(u))*(du⋅v))
  dhdrag(u,h,dh,v) = -Cd/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*(u⋅v)*dh
  dhconvₕ(u,dh,∇u,∇dh,w) = ((u⋅∇dh) + (dh)*tr(∇u))*w
  duabsᵤ(u,du) = 1.0 / (u⋅u + 1.0e-8).^(1/2) * (u⋅du)

  # Stabilization
  c₁ = 12.0; c₂ = 2.0; c₃ = 1.0
  Rₕ(u,h,hₜ,∇u,∇h) = hₜ + u⋅∇h + (h+h₀⬇)*(tr(∇u))
  Rᵤ(u,h,uₜ,∇u,∇h) = uₜ + u⋅∇u + g*∇h + Cd/(h+h₀⬇)*(absᵤ(u))*u
  Lᵤᵃ(u,∇v,∇w) = - u⋅∇v #- g*∇w
  Lₕᵃ(u,h,∇v,∇w) = - u⋅∇w - (h+h₀⬇)*(tr(∇v))
  τᵤinv(a,h,Δx₀) = 1.0/Δt + (c₁*ν / (Δx₀*Δx₀)) + (c₂*absᵤ(a) / Δx₀) + (c₃*Cd*g*absᵤ(a) / (h+1.0e-8))
  τᵤ(a,h,Δx₀) = 1.0 / τᵤinv(a,h,Δx₀)
  τₕ(a,h,Δx₀) = (Δx₀^2)/(c₁*τᵤ(a,h,Δx₀))
  stabₕ(u,h,hₜ,∇u,∇h,∇v,∇w,Δx₀) = (τₕ(u,h,Δx₀))*((Rₕ(u,h,hₜ,∇u,∇h))*Lₕᵃ(u,h,∇v,∇w))
  stabᵤ(u,h,uₜ,∇u,∇h,∇v,∇w,Δx₀) = (τᵤ(u,h,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w))
  dustabᵤ(u,du,h,uₜ,∇u,∇du,∇h,∇v,∇w,Δx₀) =
    (duτᵤ(u,du,h,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
    (τᵤ(u,h,Δx₀))*((duRᵤ(u,du,h,uₜ,∇u,∇du,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
    (τᵤ(u,h,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅duLᵤᵃ(du,∇v,∇w))
  dhstabᵤ(u,h,dh,uₜ,∇u,∇h,∇dh,∇v,∇w,Δx₀) =
    (dhτᵤ(u,h,dh,Δx₀))*((Rᵤ(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w)) +
    (τᵤ(u,h,Δx₀))*((dhRᵤ(u,h,dh,uₜ,∇u,∇dh))⋅Lᵤᵃ(u,∇v,∇w))
  duₜstabᵤ(u,h,duₜ,∇u,∇h,∇v,∇w,Δx₀) = (τᵤ(u,h,Δx₀))*((duₜRᵤ(duₜ))⋅Lᵤᵃ(u,∇v,∇w))
  duτᵤ(a,da,h,Δx₀) = -1.0 / (τᵤinv(a,h,Δx₀)*τᵤinv(a,h,Δx₀)) * duτᵤinv(a,da,h,Δx₀)
  duτᵤinv(a,da,h,Δx₀) = (c₂*duabsᵤ(a,da) / Δx₀) + (c₃*Cd*g*duabsᵤ(a,da) / (h+1.0e-8))
  dhτᵤ(a,h,dh,Δx₀) = -1.0 / (τᵤinv(a,h,Δx₀)*τᵤinv(a,h,Δx₀)) * dhτᵤinv(a,h,dh,Δx₀)
  dhτᵤinv(a,h,dh,Δx₀) = -1.0*(c₃*Cd*g*absᵤ(a) / ((h+1.0e-8)*(h+1.0e-8))) * dh
  duRᵤ(u,du,h,uₜ,∇u,∇du,∇h) = u⋅∇du + du⋅∇u + Cd/(h+h₀⬇)*(duabsᵤ(u,du))*u + Cd/(h+h₀⬇)*(absᵤ(u))*du
  duLᵤᵃ(du,∇v,∇w) = - du⋅∇v
  dhRᵤ(u,h,dh,uₜ,∇u,∇dh) = g*∇dh - Cd/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*u *dh
  duₜRᵤ(duₜ) = duₜ


  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = dΩ.trian
  # Ω = get_triangulation(dΩ.quad)
  # Δx₀ = lazy_map(dx->dx^(1/D),get_cell_measure(Ω))
  Δx₀ = CellField(_get_cell_size(Ω),Ω)

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        (strs∘(∇(u),∇(v))) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) -
                        (stabᵤ(u,h,∂t(u),∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ#-
                        # (stabₕ(u,h,∂t(h),∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))
  jac(t,(u,h),(du,dh),(v,w)) =
    ∫( (convᵤ(du,∇(u),v))  )dΩ +
    ∫( (convᵤ(u,∇(du),v)) )dΩ +
    ∫( (strs(∇(du),∇(v))) )dΩ +
    ∫( (dudrag(u,du,h,v)) )dΩ +
    ∫( (dhdrag(u,h,dh,v)) )dΩ +
    ∫( (grad(∇(dh),v)) )dΩ +
    ∫( (convₕ(du,h,∇(du),∇(h),w)) )dΩ +
    ∫( (dhconvₕ(u,dh,∇(u),∇(dh),w)) )dΩ -
    ∫( (dustabᵤ(u,du,h,∂t(u),∇(u),∇(du),∇(h),∇(v),∇(w),Δx₀)) )dΩ -
    ∫( (dhstabᵤ(u,h,dh,∂t(u),∇(u),∇(h),∇(dh),∇(v),∇(w),Δx₀)) )dΩ
  jac_t(t,(u,h),(duₜ,dhₜ),(v,w)) =
    m(t,(duₜ,dhₜ),(v,w)) -
    ∫( (duₜstabᵤ(u,h,duₜ,∇(u),∇(h),∇(v),∇(w),Δx₀)) )dΩ

  return m,a,res,jac,jac_t

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
  νₜ(εᵤ,Δx₀) = (cₛ^2*(Δx₀*Δx₀))*(√(2*(εᵤ⊙εᵤ)+1.0e-8 ))
  absᵤ(u) = (u⋅u + 1.0e-8).^(1/2)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  # strs(∇u,∇v,εᵤ,Δx₀) = ( (ν+νₜ(εᵤ,Δx₀))*(∇u+∇u') - 2/3*(ν+νₜ(εᵤ,Δx₀))*tr(∇u)*I) ⊙ ∇v
  strs(ν,∇u,∇v) = ( (ν)*(∇u+∇u') - 2/3*(ν)*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd^2/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(v⋅∇h)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w

  # Derivatives
  dabsᵤ(u,du) = 1/absᵤ(u)*(du⋅u)
  dudrag(u,du,h,v) = Cd^2/(h+h₀⬇)*((dabsᵤ(u,du))*(u⋅v)+(absᵤ(u))*(du⋅v))
  dhdrag(u,h,dh,v) = -Cd^2/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*(u⋅v)*dh
  dhconvₕ(u,dh,∇u,∇dh,w) = ((u⋅∇dh) + (dh)*tr(∇u))*w
  νₜ(εᵤ,εdu,Δx₀) = (cₛ^2*Δx₀*Δx₀)/((2*(εᵤ⊙εᵤ)+1.0e-8 ).^(1/2 ))*(2*(εdu⊙εᵤ))

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = get_triangulation(dΩ.quad)
  Δx₀ = lazy_map(dx->dx^(1/D),get_cell_measure(Ω))

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        # (strs∘(∇(u),∇(v),ε(u),Δx₀)) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ +
                     ∫( (strs(ν,∇(u),∇(v))) )dΩ +
                     ∫( (strs(νₜ∘(ε(u),Δx₀),∇(u),∇(v))) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))
  jac(t,(u,h),(du,dh),(v,w)) =
    ∫( (convᵤ(du,∇(u),v))  )dΩ +
    ∫( (convᵤ(u,∇(du),v)) )dΩ +
    ∫( (strs(ν,∇(du),∇(v))) )dΩ +
    ∫( (strs(νₜ∘(ε(u),Δx₀),∇(du),∇(v))) )dΩ +
    ∫( (strs(νₜ∘(ε(u),ε(du),Δx₀),∇(u),∇(v))) )dΩ +
    ∫( (dudrag(u,du,h,v)) )dΩ +
    ∫( (dhdrag(u,h,dh,v)) )dΩ +
    ∫( (grad(∇(dh),v)) )dΩ +
    ∫( (convₕ(du,h,∇(du),∇(h),w)) )dΩ +
    ∫( (dhconvₕ(u,dh,∇(u),∇(dh),w)) )dΩ
  jac_t(t,(uₜ,hₜ),(duₜ,dhₜ),(v,w)) = m(t,(duₜ,dhₜ),(v,w))

  return m,a,res,(jac,jac_t)

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
  νₜ(εᵤ,Δx₀) = (cₛ^2*(Δx₀*Δx₀))*(√(2*(εᵤ⊙εᵤ)+1.0e-8 ))
  absᵤ(u) = (u⋅u + 1.0e-8).^(1/2)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  # strs(∇u,∇v,εᵤ,Δx₀) = ( (ν+νₜ(εᵤ,Δx₀))*(∇u+∇u') - 2/3*(ν+νₜ(εᵤ,Δx₀))*tr(∇u)*I) ⊙ ∇v
  strs(ν,∇u,∇v) = ( (ν)*(∇u+∇u') - 2/3*(ν)*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd^2/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(v⋅∇h)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w

  # Derivatives
  dabsᵤ(u,du) = 1/absᵤ(u)*(du⋅u)
  dudrag(u,du,h,v) = Cd^2/(h+h₀⬇)*((dabsᵤ(u,du))*(u⋅v)+(absᵤ(u))*(du⋅v))
  dhdrag(u,h,dh,v) = -Cd^2/((h+h₀⬇)*(h+h₀⬇))*(absᵤ(u))*(u⋅v)*dh
  dhconvₕ(u,dh,∇u,∇dh,w) = ((u⋅∇dh) + (dh)*tr(∇u))*w
  νₜ(εᵤ,εdu,Δx₀) = (cₛ^2*Δx₀*Δx₀)/((2*(εᵤ⊙εᵤ)+1.0e-8 ).^(1/2 ))*(2*(εdu⊙εᵤ))

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = dΩ.trian
  Δx₀ = CellField(_get_cell_size(Ω),Ω)

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        # (strs∘(∇(u),∇(v),ε(u),Δx₀)) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ +
                     ∫( (strs(ν,∇(u),∇(v))) )dΩ +
                     ∫( (strs(νₜ∘(ε(u),Δx₀),∇(u),∇(v))) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))
  jac(t,(u,h),(du,dh),(v,w)) =
    ∫( (convᵤ(du,∇(u),v))  )dΩ +
    ∫( (convᵤ(u,∇(du),v)) )dΩ +
    ∫( (strs(ν,∇(du),∇(v))) )dΩ +
    ∫( (strs(νₜ∘(ε(u),Δx₀),∇(du),∇(v))) )dΩ +
    ∫( (strs(νₜ∘(ε(u),ε(du),Δx₀),∇(u),∇(v))) )dΩ +
    ∫( (dudrag(u,du,h,v)) )dΩ +
    ∫( (dhdrag(u,h,dh,v)) )dΩ +
    ∫( (grad(∇(dh),v)) )dΩ +
    ∫( (convₕ(du,h,∇(du),∇(h),w)) )dΩ +
    ∫( (dhconvₕ(u,dh,∇(u),∇(dh),w)) )dΩ
  jac_t(t,(uₜ,hₜ),(duₜ,dhₜ),(v,w)) = m(t,(duₜ,dhₜ),(v,w))

  return m,a,res,(jac,jac_t)

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

# Auxiliar functions
function _get_cell_size(t::Triangulation)
  meas = get_cell_measure(t)
  d = num_dims(t)
  map(m->m^(1/d),meas)
end

function _get_cell_size(t::GridapDistributed.DistributedTriangulation)
  map(_get_cell_size,local_views(t))
end

"""
    get_forms(::Val{:Galerkin},params::physics_params)

Get the FE operator forms for the standard Galerkin method.
"""
function get_forms(measures,normals,D,::Val{:conservative_Galerkin},
  physics_params::physics_params,
  ode_solver_params::ODE_solver_params)

  @unpack ν,Cd,g,h₀⬇ = physics_params
  if D==1
    @error "1D not implemented"
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  H(h) = h+h₀⬇
  𝒜₁(u) = TensorValue(
    0, 1.0, 0,
    g*u[1]-(u[2]/u[1])^2, 2*u[2]/u[1], 0,
    -u[2]*u[3]/(u[1]*u[1]), u[3]/u[1], u[2]/u[1]
  )
  𝒜₂(u) = TensorValue(
    0, 0, 1.0,
    -u[2]*u[3]/(u[1]*u[1]), u[3]/u[1], u[2]/u[1],
    g*u[1]-(u[3]/u[1])^2, 0, 2*u[3]/u[1]
  )
  𝒦₁₁(u) = TensorValue(
    0, 0, 0,
    -2ν*u[2]/u[1], 2ν, 0,
    -ν*u[3]/u[1], 0, ν
  )
  𝒦₁₂(u) = TensorValue(
    0, 0, 0,
    0, 0, 0,
    -ν*u[2]/u[1], ν, 0
  )
  𝒦₂₁(u) = TensorValue(
    0, 0, 0,
    -ν*u[3]/u[1], 0, ν,
    0, 0, 0
  )
  𝒦₂₂(u) = TensorValue(
    0, 0, 0,
    -ν*u[2]/u[1], ν, 0,
    -2ν*u[3]/u[1], 0, 2ν
  )
  𝒮(u) = VectorValue(
    0.0,
    0.0,#- g*Cd^2*_abs(u)*u[2],#/((u[1] + 1.0e-8)^(1/3)),
    0.0#- g*Cd^2*_abs(u)*u[3]#/((u[1] + 1.0e-8)^(1/3))
  )
  ℋ(u,n) = VectorValue(
    0,
    0.5*g*u[1]^2*n[1],
    0.5*g*u[1]^2*n[2]
  )
  R(u,∇₁u,∇₂u,εu) = VectorValue(
    0,
    -2ν*(∇₁u[1]*εu[1,1]+∇₂u[1]*εu[1,2]),# + g*Cd*_abs(u)*u[2]/((u[1] + 1.0e-8)^(1/3)),
    -2ν*(∇₁u[1]*εu[2,1]+∇₂u[1]*εu[2,2])# + g*Cd*_abs(u)*u[3]/((u[1] + 1.0e-8)^(1/3))
  )
  ℛ(u) = ∂t(u) - (𝒮∘u) + ℒ(u) #- (R∘(u,∇₁(u),∇₂(u),εᵤ(u))) # only true for 1st order
  # ℛ(u) = 𝒵(u) - (R∘(u,∇₁(u),∇₂(u),εᵤ(u)))
  ℒ(u) = (𝒜₁∘u)⋅(∇₁(u)) + (𝒜₂∘u)⋅∇₂(u)
  ℒᵃ(w,u) = ∇₁(w)⋅(𝒜₁∘u) + ∇₂(w)⋅(𝒜₂∘u)

  ∇₁(u) = VectorValue(1.0,0.0)⋅∇(u)
  ∇₂(u) = VectorValue(0.0,1.0)⋅∇(u)
  ∇ᵤ(u) = ∇(u)⋅TensorValue{3,2}(0.0,0.0,1.0,0.0,0.0,1.0)
  εᵤ(u) = 1/2*(∇ᵤ(u) + ∇ᵤ(u)')

  dΩ,dΓwall, = measures
  nwall, = normals
  Ω = get_triangulation(dΩ.quad)
  h = lazy_map(dx->dx^(1/D),get_cell_measure(Ω))
  _abs(u) = √(u[2]^2+u[3]^2 + 1.0e-8)
  ∇h(∇u) = ∇u⋅VectorValue(1.0,0.0,0.0)
  _absh(∇h) = √(∇h[1]^2+∇h[2]^2 + 1.0e-8)
  τshoc(u,∇u,h) = h/(2*_abs(u))*(_absh(∇h(∇u))*h/(h₀⬇))
  νshoc(u,∇u,h) = τshoc(u,∇u,h)*_abs(u)^2

  @unpack Δt = ode_solver_params
  cτ = 0.5
  # τ = cτ*Δt/2
  τ(u,h) = 1/(2/(cτ*Δt) + 12*ν/h^2 + 2*_abs(u)/h)

  # Residual form
  dΩ,dΓwall, = measures
  nwall, = normals
  # m(t,uₜ,w) = ∫( uₜ⋅w )dΩ
  # res(t,u,w) = ∫( ℛ(u)⋅w +
  #   ((𝒦₁₁∘u)⋅(∇₁(u)) + (𝒦₁₂∘u)⋅∇₂(u))⊙(∇₁(w)) +
  #   ((𝒦₂₁∘u)⋅(∇₁(u)) + (𝒦₂₂∘u)⋅∇₂(u))⊙(∇₂(w)) +
  #   (τ∘(u,h))*((∇₁(w)⋅(𝒜₁∘u) + ∇₂(w)⋅(𝒜₂∘u))⋅ℛ(u)) )dΩ#+
    # (νshoc∘(u,∇(u),h))*(∇₁(u)⋅∇₁(w) + ∇₂(u)⋅∇₂(w)) )dΩ
  res(t,u,w) = ∫( w⋅ℛ(u) + ℒᵃ(w,u)⋅((τ∘(u,h))*(ℛ(u))) + (νshoc∘(u,∇(u),h))*(∇₁(u)⋅∇₁(w) + ∇₂(u)⋅∇₂(w)) )dΩ #+
  # res(t,u,w) = ∫( w⋅ℛ(u) + ℒᵃ(w,u)⋅((τ)*(ℛ(u))) )dΩ#+ (νshoc∘(u,∇(u),h))*(∇₁(u)⋅∇₁(w) + ∇₂(u)⋅∇₂(w)) )dΩ #+
    # ∫( w⋅(ℋ∘(u,nwall)) )dΓwall

    # ((𝒦₁₁∘u)⋅(∇₁(u)) + (𝒦₁₂∘u)⋅∇₂(u))⊙(∇₁(w)) +
    # ((𝒦₂₁∘u)⋅(∇₁(u)) + (𝒦₂₂∘u)⋅∇₂(u))⊙(∇₂(w)) +
    # (τ∘(u,h))*((∇₁(w)⋅(𝒜₁∘u) + ∇₂(w)⋅(𝒜₂∘u))⋅ℛ(u)) )dΩ#+
  # (νshoc∘(u,∇(u),h))*(∇₁(u)⋅∇₁(w) + ∇₂(u)⋅∇₂(w)) )dΩ
  # res(t,(u,h),(v,w)) = m(t,∂t(u),v) + a(t,u,w)

  return nothing,nothing,res

end

# FE operator
function get_FEOperator(forms,X,Y,::Union{Val{:Galerkin}})#,Val{:Smagorinsky}})
  m,a,res = forms
  return TransientSemilinearFEOperator(m,a,X,Y)
end
function get_FEOperator(forms,X,Y,::Union{Val{:ASGS},Val{:conservative_Galerkin}})
  # _,_,res,jacs = forms
  _,_,res = forms
  # return TransientFEOperator(res,jacs,X,Y)
  return TransientFEOperator(res,X,Y)
end
function get_FEOperator(forms,X,Y,::Val{:Smagorinsky})
  m,a,res,jacs = forms
  return TransientSemilinearFEOperator(m,a,jacs,X,Y)
end
