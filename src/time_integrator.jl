abstract type ODE_solver_params end

@with_kw struct theta_method_params <: ODE_solver_params
  type::Symbol = :θ_method
  θ::Float64 = 1.0
  Δt::Float64 = 0.25
  T::Float64 = 0.25
end

@with_kw struct Generalized_α_params <: ODE_solver_params
  type::Symbol = :Generalized_α
  ρ∞::Float64 = 1.0
  Δt::Float64 = 0.25
  T::Float64 = 0.25
end

@with_kw struct EXRK_params <: ODE_solver_params
  type::Symbol = :EXRK_SSP_3_3
  Δt::Float64 = 0.25
  T::Float64 = 0.25
end

ODE_solver_params() = theta_method_params()
function ODE_solver_params(type::Symbol;kwargs...)
  if type == :θ_method
    return theta_method_params(;kwargs...)
  elseif type == :Generalized_α
    return Generalized_α_params(;kwargs...)
  elseif type == :EXRK_SSP_3_3
    return EXRK_params(;kwargs...)
  else
    error("ODE solver type not recognized")
  end
end

function get_ode_solver(nls,params::theta_method_params)
  @unpack type,T,Δt,θ = params
  return ThetaMethod(nls,Δt,θ)
end

function get_ode_solver(nls,params::Generalized_α_params)
  @unpack type,T,Δt,ρ∞ = params
  return GeneralizedAlpha1(nls,Δt,ρ∞)
end

function get_ode_solver(nls,params::EXRK_params)
  @unpack type,T,Δt = params
  return RungeKutta(nls.ls,nls.ls,Δt,type)
end

function get_solution(odes,op,x₀,xdot₀,params::Union{theta_method_params,EXRK_params})
  @unpack T = params
  return solve(odes,op,0.0,T,x₀)
end

function get_solution(odes,op,x₀,xdot₀,params::Generalized_α_params)
  @unpack T = params
  return solve(odes,op,0.0,T,(x₀,xdot₀))
end
