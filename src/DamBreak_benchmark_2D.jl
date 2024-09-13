"""
    struct DamBreak_benchmark_2D_params

A structure to store the parameters of the 2D Dam Break benchmark. The parameters are:
  - `mesh_file::String`: Mesh file
  - `x₀::Float64`: Initial position of the dam
  - `order::Int64`: Order of the velocity FE space
  - `formulation::Symbol`: Formulation
  - `verbose::Bool`: Verbosity
  - `vtk_output::Bool`: Output VTK files
  - `vtk_folder::String`: Folder name to store VTK files
  - `physics_params::physics_params`: Physical parameters
  - `ode_solver_params::ODE_solver_params`: ODE solver parameters
"""
@with_kw struct DamBreak_benchmark_2D_params
  mesh_file::String = datadir("meshes","DamBreak_benchmark_2D_coarse.msh")  # Mesh file
  x₀::Float64 = 100.0  # Initial position of the dam
  order::Int64 = 2      # Order of the velocity FE space
  formulation::Symbol = :Galerkin # Formulation
  verbose::Bool = false # Verbosity
  vtk_output::Bool = false # Output VTK files
  vtk_folder::String = "test"  # Folder name to store VTK files (inside `data/sims`)
  physics_params::physics_params = physics_params() # Physical parameters
  ode_solver_params::ODE_solver_params = ODE_solver_params() # ODE solver parameters
end

"""
    main(params::DamBreak_benchmark_2D_params)

Main function to run the 2D Dam Break benchmark.
"""
function main(params::DamBreak_benchmark_2D_params)

  # Define the domain and Triangulations
  @unpack mesh_file, verbose = params
  model = GmshDiscreteModel(mesh_file)
  Ω = Interior(model)
  Γ = Boundary(model,tags="wall")

  # Define boundary conditions
  u₀(x,t) = VectorValue(0.0,0.0)
  u₀(t::Real) = x->u₀(x,t)

  # Define initial conditions
  @unpack x₀, physics_params = params
  @unpack h₀⬆, h₀⬇ = physics_params
  h₀(x) = x[1] < x₀ ? (h₀⬆-h₀⬇) : 0.0

  # Define spaces
  @unpack order, formulation = params
  X,Y = get_FESpaces(Ω,order,[],[],[],Val(formulation))

  # Integration Measure
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)
  measures = (dΩ,dΓ)

  # Normals
  nΓ = get_normal_vector(Γ)
  normals = (nΓ,)

  # Weak form
  @unpack ode_solver_params = params
  m,a,res = get_forms(measures,normals,2,Val(formulation), physics_params, ode_solver_params)
  op = TransientSemilinearFEOperator(m,a,X,Y)

  # Solver
  ls = LUSolver()
  nls = NLSolver(ls,show_trace=verbose,iterations=10,method=:newton)
  odes = get_ode_solver(nls,params.ode_solver_params)

  # Initial solution
  xₕ₀ = interpolate_everywhere([u₀(0),h₀],X)
  xdotₕ₀ = interpolate_everywhere([VectorValue(0.0,0.0),0.0],X)
  @unpack vtk_output, vtk_folder = params
  vtk_output && writevtk(Ω,datadir("sims","test","sol_DB2D_0.vtu"),cellfields=["u"=>xₕ₀[1],"h"=>xₕ₀[2]])

  # Solution
  xₕₜ = get_solution(odes,op,xₕ₀,xdotₕ₀,params.ode_solver_params)

  # Iterate over time
  @unpack T = params.ode_solver_params
  createpvd(datadir("sims",vtk_folder,"sol_DB2D")) do pvd
    for (t,(uₕ,hₕ)) in xₕₜ
      println("Time: $t / $T")
      vtk_output && (pvd[t] = createvtk(Ω,datadir("sims",vtk_folder,"sol_DB2D_$(t).vtu"),cellfields=["u"=>uₕ,"h"=>hₕ],order=order))
    end
  end

  return nothing

end
