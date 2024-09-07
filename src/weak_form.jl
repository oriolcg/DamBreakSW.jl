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

# Residual functions
"""
    get_residual_form(::Val{:Galerkin},params::physics_params)

Get the residual form for the standard Galerkin method.
"""
function get_residual_form(dΩ,D,::Val{:Galerkin},params::physics_params)

  @unpack ν,Cd,g,h₀⬇ = params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  absᵤ(u) = √(u⋅u + 1.0e-14)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  strs(∇u,∇v) = ( ν*(∇u+∇u') - 2/3*ν*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(v⋅∇h)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        (strs∘(∇(u),∇(v))) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))

  return res

end


"""
    get_residual_form(::Val{:ASGS},params::physics_params)

Get the residual form for the stabilized formulation using Algebraic Subgrid Scales (ASGS).
"""
function get_residual_form(dΩ,D,::Val{:ASGS},params::physics_params)

  @unpack ν,Cd,g,h₀⬇ = params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  absᵤ(u) = (u⋅u + 1.0e-14).^(1/2)
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
  τᵤ(a,h,Δxₒ) = 1.0 / ((c₁*ν / (Δxₒ^2)) + (c₂*absᵤ(a) ./ Δxₒ) + (c₃*Cd*g*absᵤ(a) / (h+1e-14)))
  τₕ(a,h,Δxₒ) = (Δxₒ^2)/(c₁*τᵤ(a,h,Δxₒ))
  stabₕ(u,h,hₜ,∇u,∇h,∇v,∇w,Δxₒ) = (τₕ∘(u,h,Δxₒ))*((Rₕ∘(u,h,hₜ,∇u,∇h))*Lₕᵃ(u,h,∇v,∇w))
  stabᵤ(u,h,uₜ,∇u,∇h,∇v,∇w,Δxₒ) = (τᵤ∘(u,h,Δxₒ))*((Rᵤ∘(u,h,uₜ,∇u,∇h))⋅Lᵤᵃ(u,∇v,∇w))

  Ω = get_triangulation(dΩ.quad)
  Δxₒ = lazy_map(dx->dx^(1/D),get_cell_measure(Ω))

  # Residual form
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        (strs∘(∇(u),∇(v))) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) -
                        (stabᵤ(u,h,∂t(u),∇(u),∇(h),∇(v),∇(w),Δxₒ)) -
                        (stabₕ(u,h,∂t(h),∇(u),∇(h),∇(v),∇(w),Δxₒ)) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))

  return res

end


"""
    get_residual_form(::Val{:Smagorinsky},params::physics_params)

Get the residual form for the standard Galerkin method.
"""
function get_residual_form(dΩ,D,::Val{:Smagorinsky},params::physics_params)

  @unpack ν,Cd,g,h₀⬇ = params
  if D==1
    I = TensorValue(1.0)
  elseif D==2
    I = TensorValue(1.0,0.0,0.0,1.0)
  end

  # Auxiliar functions
  cₛ = 0.2
  νₜ(εᵤ,Δx₀) = cₛ*Δx₀^2*(√(2*(εᵤ⊙εᵤ)+1.0e-14))
  absᵤ(u) = √(u⋅u + 1.0e-14)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  strs(∇u,∇v,εᵤ,Δx₀) = ( (ν+νₜ(εᵤ,Δx₀))*(∇u+∇u') - 2/3*(ν+νₜ(εᵤ,Δx₀))*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(v⋅∇h)
  convₕ(u,h,∇u,∇h,w) = ((u⋅∇h) + (h+h₀⬇)*tr(∇u))*w

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

  return res

end
