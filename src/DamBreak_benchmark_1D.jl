"""
    struct DamBreak_benchmark_1D_params

A structure to store the parameters of the 1D Dam Break benchmark. The parameters are:
  - `L::Float64`: Length of the domain
  - `x₀::Float64`: Initial position of the dam
  - `nₓ::Int64`: Number of cells
  - `order::Int64`: Order of the velocity FE space
  - `formulation::Symbol`: Formulation
  - `verbose::Bool`: Verbosity
  - `vtk_output::Bool`: Output VTK files
  - `vtk_folder::String`: Folder name to store VTK files
  - `physics_params::physics_params`: Physical parameters
  - `ode_solver_params::ODE_solver_params`: ODE solver parameters
"""
@with_kw struct DamBreak_benchmark_1D_params
  L::Float64 = 2000.0   # Length of the domain
  x₀::Float64 = 1000.0  # Initial position of the dam
  nₓ::Int64 = 4         # Number of cells
  order::Int64 = 2      # Order of the velocity FE space
  formulation::Symbol = :Galerkin # Formulation
  verbose::Bool = false # Verbosity
  vtk_output::Bool = false # Output VTK files
  vtk_folder::String = "test"  # Folder name to store VTK files (inside `data/sims`)
  physics_params::physics_params = physics_params() # Physical parameters
  ode_solver_params::ODE_solver_params = ODE_solver_params() # ODE solver parameters
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
  Γ = Boundary(model)
  Λ = Skeleton(model)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"outlet",[2])
  Γout = Boundary(model,tags="outlet")

  # Define boundary conditions
  u₀(x,t) = VectorValue(0.0)
  u₀(t::Real) = x->u₀(x,t)

  # Define initial conditions
  @unpack x₀, physics_params = params
  @unpack h₀⬆, h₀⬇ = physics_params
  h₀(x) = x[1] < x₀ ? (h₀⬆-h₀⬇) : 0.0

  # Define spaces
  @unpack order, formulation = params
  X,Y = get_FESpaces(Ω,1,order,["boundary"],[(true,)],[u₀],Val(formulation))

  # Integration Measure
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)
  dΛ = Measure(Λ,2*order)
  dΓout = Measure(Γout,2*order)
  measures = (dΩ,dΓ,dΛ,dΓout)

  # Normals
  normals = (get_normal_vector(Γ),get_normal_vector(Λ),get_normal_vector(Γout))

  # Weak form
  @unpack ode_solver_params = params
  forms = get_forms(measures,normals,1,Val(formulation), physics_params, ode_solver_params)
  op = get_FEOperator(forms,X,Y,Val(formulation))

  # Solver
  ls = LUSolver()
  nls = NLSolver(ls,show_trace=verbose,iterations=10,method=:newton)
  odes = get_ode_solver(nls,ode_solver_params)

  # Initial solution
  xₕ₀ = interpolate_everywhere([u₀(0),h₀],X)
  xdotₕ₀ = interpolate_everywhere([VectorValue(0.0),0.0],X)
  @unpack vtk_output, vtk_folder = params
  vtk_output && writevtk(Ω,datadir("sims","test","sol_DB1D_0.vtu"),cellfields=["u"=>xₕ₀[1],"h"=>xₕ₀[2]])

  # Solution
  xₕₜ = get_solution(odes,op,xₕ₀,xdotₕ₀,ode_solver_params)

  # Iterate over time
  @unpack T = params.ode_solver_params
  createpvd(datadir("sims",vtk_folder,"sol_DB1D")) do pvd
    for (t,(uₕ,hₕ)) in xₕₜ
      println("Time: $t / $T")
      vtk_output && (pvd[t] = createvtk(Ω,datadir("sims",vtk_folder,"sol_DB1D_$(t).vtu"),cellfields=["u"=>uₕ,"h"=>hₕ]))
    end
  end

  return nothing

end
