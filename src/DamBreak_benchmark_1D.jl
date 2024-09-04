"""
    struct DamBreak_benchmark_1D_params

A structure to store the parameters of the 1D Dam Break benchmark. The parameters are:
  - `L::Float64`: Length of the domain
  - `x₀::Float64`: Initial position of the dam
  - `h₀⬆::Float64`: Initial height of the water upstream
  - `h₀⬇::Float64`: Initial height of the water downstream
  - `nₓ::Int64`: Number of cells
  - `order::Int64`: Order of the velocity FE space
  - `T::Float64`: Final time
  - `Δt::Float64`: Time step
  - `ν::Float64`: Viscosity
  - `Cd::Float64`: Drag coefficient
  - `verbose::Bool`: Verbosity
  - `vtk_output::Bool`: Output VTK files
  - `vtk_folder::String`: Folder name to store VTK files
"""
@with_kw struct DamBreak_benchmark_1D_params
  L::Float64 = 2000.0   # Length of the domain
  x₀::Float64 = 1000.0  # Initial position of the dam
  h₀⬆::Float64 = 10.0   # Initial height of the water upstream
  h₀⬇::Float64 = 5.0    # Initial height of the water downstream
  nₓ::Int64 = 4         # Number of cells
  order::Int64 = 2      # Order of the velocity FE space
  T::Float64 = 0.25     # Final time
  Δt::Float64 = 0.25    # Time step
  θ::Float64 = 1.0      # θ-scheme
  ν::Float64 = 1.0      # Viscosity
  Cd::Float64 = 0.0     # Drag coefficient
  g::Float64 = 9.8      # Gravity
  verbose::Bool = false # Verbosity
  vtk_output::Bool = false # Output VTK files
  vtk_folder::String = "test"  # Folder name to store VTK files (inside `data/sims`)
end

"""
    main(params::DamBreak_benchmark_1D_params)

Main function to run the 1D Dam Break benchmark.
"""
function main(params::DamBreak_benchmark_1D_params)

  # Define the domain and Triangulations
  @unpack L, nₓ, verbose = params
  model = CartesianDiscreteModel((0,L), (nₓ,))
  Ω = Interior(model)

  # Define boundary conditions
  u₀(x,t) = VectorValue(0.0)
  u₀(t::Real) = x->u₀(x,t)

  # Define initial conditions
  @unpack x₀, h₀⬆, h₀⬇ = params
  h₀(x) = x[1] < x₀ ? (h₀⬆-h₀⬇) : 0.0

  # Define spaces
  @unpack order = params
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
  refFEₕ = ReferenceFE(lagrangian,Float64,order-1)
  Vᵤ = TestFESpace(Ω,refFEᵤ,dirichlet_tags="boundary")
  Vₕ = TestFESpace(Ω,refFEₕ;conformity=:H1)
  Uᵤ = TransientTrialFESpace(Vᵤ,u₀)
  Uₕ = TransientTrialFESpace(Vₕ)
  Y = MultiFieldFESpace([Vᵤ,Vₕ])
  X = TransientMultiFieldFESpace([Uᵤ,Uₕ])

  # Integration Measure
  dΩ = Measure(Ω,2*order)

  # Weak form
  @unpack ν,Cd,g = params
  I = TensorValue(1.0)
  absᵤ(u) = √(u⋅u + 1.0e-14)
  convᵤ(a,∇u,v) = (a⋅∇u)⋅v
  strs(∇u,∇v) = ( ν*(∇u+∇u') - 2/3*ν*tr(∇u)*I) ⊙ ∇v
  drag(u,h,v) = Cd/(h+h₀⬇)*(absᵤ(u))*(u⋅v)
  grad(∇h,v) = g*(∇h⋅v)
  convₕ(u,h,∇u,∇h,w) = ((∇h⋅u) + (h+h₀⬇)*tr(∇u))*w
  m(t,(uₜ,hₜ),(v,w)) = ∫(uₜ⋅v + hₜ*w)dΩ
  a(t,(u,h),(v,w)) = ∫( (convᵤ∘(u,∇(u),v)) +
                        (strs∘(∇(u),∇(v))) +
                        (drag∘(u,h,v)) +
                        (grad∘(∇(h),v)) +
                        (convₕ∘(u,h,∇(u),∇(h),w)) )dΩ
  res(t,(u,h),(v,w)) = m(t,(∂t(u),∂t(h)),(v,w)) + a(t,(u,h),(v,w))

  op = TransientFEOperator(res,X,Y)

  # Solver
  @unpack T, Δt,θ = params
  ls = LUSolver()
  nls = NLSolver(ls,show_trace=verbose,iterations=10,method=:newton)
  odes = ThetaMethod(nls,Δt,θ)

  # Initial solution
  xₕ₀ = interpolate_everywhere([u₀(0),h₀],X)
  @unpack vtk_output, vtk_folder = params
  vtk_output && writevtk(Ω,datadir("sims","test","sol0.vtu"),cellfields=["u"=>xₕ₀[1],"h"=>xₕ₀[2]])

  # Solution
  xₕₜ = solve(odes,op,0.0,T,xₕ₀)

  # Iterate over time
  createpvd(datadir("sims",vtk_folder,"sol")) do pvd
    for (t,(uₕ,hₕ)) in xₕₜ
      println("Time: $t / $T")
      vtk_output && (pvd[t] = createvtk(Ω,datadir("sims",vtk_folder,"sol_$(t).vtu"),cellfields=["u"=>uₕ,"h"=>hₕ]))
    end
  end

  return nothing

end
