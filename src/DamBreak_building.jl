"""
    struct DamBreak_building_params

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
@with_kw struct DamBreak_building_params
  mesh_file::String = datadir("meshes","DamBreak_building_0.0_coarse.msh")  # Mesh file
  x₀::Float64 = 0.0  # Initial position of the dam
  order::Int64 = 2      # Order of the velocity FE space
  formulation::Symbol = :Galerkin # Formulation
  verbose::Bool = false # Verbosity
  vtk_output::Bool = false # Output VTK files
  vtk_folder::String = "test"  # Folder name to store VTK files (inside `data/sims`)
  Δtout::Float64 = 0.05 # Time step for VTK output
  physics_params::physics_params = physics_params(g=9.81,ν=1.0e-6,h₀⬆=0.63,h₀⬇=0.03,Cd=0.0127) # Physical parameters
  ode_solver_params::ODE_solver_params = ODE_solver_params() # ODE solver parameters
end

"""
    main(params::DamBreak_building_params)

Main function to run the 2D Dam Break benchmark.
"""
function main(params::DamBreak_building_params)

  # Define the domain and Triangulations
  @unpack mesh_file, verbose = params
  model = GmshDiscreteModel(mesh_file)
  Ω = Interior(model)
  Γ = Boundary(model,tags=["walls"])

  # Define boundary conditions
  u₀(x,t) = VectorValue(0.0,0.0)
  u₀(t::Real) = x->u₀(x,t)

  # Define initial conditions
  @unpack x₀, physics_params = params
  @unpack h₀⬆, h₀⬇ = physics_params
  ϵ = 0.05
  h₀(x) = (0.5 * (1 + tanh((x₀-x[1]) / ϵ)))*(h₀⬆-h₀⬇)
  # h₀(x) = x[1] < x₀ ? (h₀⬆-h₀⬇) : 0.0

  # Define spaces
  @unpack order, formulation = params
  D = 2
  dirichlet_tags, dirichlet_masks, dirichlet_values = _get_dirichlet_DB_building(Val(formulation))
  X,Y = get_FESpaces(Ω,D,order,dirichlet_tags, dirichlet_masks, dirichlet_values,Val(formulation))

  # Integration Measure
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)
  measures = (dΩ,dΓ)

  # Normals
  nΓ = get_normal_vector(Γ)
  normals = (nΓ,)

  # Weak form
  @unpack ode_solver_params = params
  forms = get_forms(measures,normals,D,Val(formulation), physics_params, ode_solver_params)
  op = get_FEOperator(forms,X,Y,Val(formulation))

  # Solver
  ls = LUSolver()
  nls = NLSolver(ls,show_trace=verbose,iterations=15,method=:newton,ftol=1.0e-6)
  odes = get_ode_solver(nls,params.ode_solver_params)

  # Initial solution
  xₕ₀, xdotₕ₀ = _get_initial_solution_DB_building(u₀,h₀,h₀⬇,X,Val(formulation))
  @unpack vtk_output, vtk_folder, Δtout = params
  vtk_output && _writevtk_DB_building(Ω,vtk_folder,xₕ₀,Val(formulation))

  # Solution
  xₕₜ = get_solution(odes,op,xₕ₀,xdotₕ₀,params.ode_solver_params)

  # Iterate over time
  @unpack T = params.ode_solver_params
  tout = 0.0
  createpvd(datadir("sims",vtk_folder,"sol_DB2D")) do pvd
    for (t,xₕ) in xₕₜ
      println("Time: $t / $T")
      if t >= tout
        vtk_output && (pvd[t] = _createvtk_DB_building(Ω,vtk_folder,order,xₕ,t,Val(formulation)))
        tout += Δtout
      end
    end
  end

  return nothing

end



"""
    main(ranks,params::DamBreak_building_params)

Main function to run the 2D dambreak problem with an obstacle.
"""
function main(ranks,params::DamBreak_building_params)

  # Define the domain and Triangulations
  @unpack mesh_file, verbose = params
  model = GmshDiscreteModel(ranks,mesh_file)
  Ω = Interior(model)
  Γ = Boundary(model,tags=["walls"])

  # Define boundary conditions
  u₀(x,t) = VectorValue(0.0,0.0)
  u₀(t::Real) = x->u₀(x,t)

  # Define initial conditions
  @unpack x₀, physics_params = params
  @unpack h₀⬆, h₀⬇ = physics_params
  ϵ = 0.05
  h₀(x) = (0.5 * (1 + tanh((x[1]-x₀) / ϵ)))*(h₀⬆-h₀⬇)#x[1] < x₀ ? (h₀⬆-h₀⬇) : 0.0

  # Define spaces
  @unpack order, formulation = params
  D = 2
  X,Y = get_FESpaces(Ω,D,order,["walls","inlet","sides"],[(true,true),(true,false),(false,true)],[u₀,u₀,u₀],Val(formulation))

  # Integration Measure
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)
  measures = (dΩ,dΓ)

  # Normals
  nΓ = get_normal_vector(Γ)
  normals = (nΓ,)

  # Weak form
  @unpack ode_solver_params = params
  m,a,res,jac,jac_t = get_forms(measures,normals,D,Val(formulation), physics_params, ode_solver_params)
  # op = TransientFEOperator(res,X,Y)
  op = TransientSemilinearFEOperator(m,a,(jac,jac_t),X,Y)

  # Solver
  # ls = LUSolver()
  # nls = NLSolver(ls,show_trace=verbose,iterations=10,method=:newton)
  nls = PETScNonlinearSolver()
  odes = get_ode_solver(nls,params.ode_solver_params)

  # Initial solution
  xₕ₀ = interpolate_everywhere([u₀(0),h₀],X)
  xdotₕ₀ = interpolate_everywhere([VectorValue(0.0,0.0),0.0],X)
  @unpack vtk_output, vtk_folder = params
  vtk_output && _writevtk_DB_building(Ω,vtk_folder,xₕ₀,Val(formulation))

  # Solution
  xₕₜ = get_solution(odes,op,xₕ₀,xdotₕ₀,params.ode_solver_params)

  # Iterate over time
  @unpack T = params.ode_solver_params
  tout = 0.0
  Δtout = 0.05
  createpvd(ranks,datadir("sims",vtk_folder,"sol_DB2D")) do pvd
    for (t,xₕ) in xₕₜ
      println("Time: $t / $T")
      if t >= tout
        vtk_output && (pvd[t] = _createvtk_DB_building(Ω,vtk_folder,order,xₕ,t,Val(formulation)))
        tout += Δtout
      end
    end
  end

  return nothing

end


# Auxiliary functions
function _get_dirichlet_DB_building(::Union{Val{:Galerkin},Val{:ASGS},Val{:Smagorinsky}})
  u₀(x,t) = VectorValue(0.0,0.0)
  u₀(t::Real) = x->u₀(x,t)
  return ["walls","inlet","sides"],[(true,true),(true,false),(false,true)],[u₀,u₀,u₀]
end
function _get_dirichlet_DB_building(::Val{:conservative_Galerkin})
  U₀(x,t) = VectorValue(0.0,0.0,0.0)
  U₀(t) = x->U₀(x,t)
  return ["walls","inlet","sides"], [(false,true,true),(false,true,false),(false,false,true)], [U₀,U₀,U₀]
end

function _get_initial_solution_DB_building(u₀,h₀,h₀⬇,X,::Union{Val{:Galerkin},Val{:ASGS},Val{:Smagorinsky}})
  return interpolate_everywhere([u₀(0),h₀],X), interpolate_everywhere([VectorValue(0.0,0.0),0.0],X)
end
function _get_initial_solution_DB_building(u₀,h₀,h₀⬇,X,::Val{:conservative_Galerkin})
  U₀(x) = VectorValue(h₀(x)+h₀⬇,u₀(x,0.0)[1],u₀(x,0.0)[2])
  return interpolate_everywhere(U₀,X), interpolate_everywhere(VectorValue(0.0,0.0,0.0),X)
end

function _writevtk_DB_building(Ω,vtk_folder,xₕ₀,::Union{Val{:Galerkin},Val{:ASGS},Val{:Smagorinsky}})
  writevtk(Ω,datadir("sims",vtk_folder,"sol_DB2D_0.vtu"),cellfields=["u"=>xₕ₀[1],"h"=>xₕ₀[2]])
end
function _writevtk_DB_building(Ω,vtk_folder,xₕ₀,::Val{:conservative_Galerkin})
  writevtk(Ω,datadir("sims",vtk_folder,"sol_DB2D_0.vtu"),cellfields=["U"=>xₕ₀])
end

function _createvtk_DB_building(Ω,vtk_folder,order,xₕ,t,::Union{Val{:Galerkin},Val{:ASGS},Val{:Smagorinsky}})
  createvtk(Ω,datadir("sims",vtk_folder,"sol_DB2D_$(t).vtu"),cellfields=["u"=>xₕ[1],"h"=>xₕ[2]],order=order)
end
function _createvtk_DB_building(Ω,vtk_folder,order,xₕ,t,::Val{:conservative_Galerkin})
  createvtk(Ω,datadir("sims",vtk_folder,"sol_DB2D_$(t).vtu"),cellfields=["U"=>xₕ],order=order)
end
