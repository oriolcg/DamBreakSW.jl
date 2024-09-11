"""
    struct Mesh_params

Parameters for the mesh generation. The following parameters are considered:
  - `x₁::Vector{2,Float64}`: Lower left corner
  - `x₂::Vector{2,Float64}`: Upper right corner
  - `x₃::Vector{2,Float64}`: Lower left corner of the building (rectangular)
  - `x₄::Vector{2,Float64}`: Upper right corner of the building (rectangular)
  - `H::Float64`: Coarse element size
  - `h::Float64`: Fine element size
  - `dxLeft::Float64`: Refinement zone distance from object (left)
  - `dxRight::Float64`: Refinement zone distance from object (right)
  - `dyTop::Float64`: Refinement zone distance from object (top)
  - `dyBottom::Float64`: Refinement zone distance from object (bottom)
  - `decay_factor_left::Float64`: Decay factor for the left refinement zone
  - `decay_factor_right::Float64`: Decay factor for the right refinement zone
  - `decay_exponent_left::Float64`: Decay exponent for the left refinement zone
  - `decay_exponent_right::Float64`: Decay exponent for the right refinement zone
"""
@with_kw struct Mesh_params
  x₁::Vector{2,Float64} = [-11.678,0.0]  # Lower left corner
  x₂::Vector{2,Float64} = [32.0,1.4]     # Upper right corner
  x₃::Vector{2,Float64} = [14.0,0.55]    # Lower left corner of the building (rectangular)
  x₄::Vector{2,Float64} = [14.3,0.85]    # Upper right corner of the building (rectangular)
  H::Float64 = 0.5                       # Coarse element size
  h::Float64 = 0.1                       # Fine element size
  dxLeft::Float64 = 0.5                  # Refinement zone distance from object (left)
  dxRight::Float64 = 6.0                 # Refinement zone distance from object (right)
  dyTop::Float64 = 0.5                   # Refinement zone distance from object (top)
  dyBottom::Float64 = 0.5                # Refinement zone distance from object (bottom)
  decay_factor_left::Float64 = 5.0       # Decay factor for the left refinement zone
  decay_factor_right::Float64 = 1.5      # Decay factor for the right refinement zone
  decay_exponent_left::Float64 = 0.8     # Decay exponent for the left refinement zone
  decay_exponent_right::Float64 = 0.8    # Decay exponent for the right refinement zone
end


function create_mesh_length(params::Mesh_params)

  # Initialize
  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  gmsh.model.add("t1")

  # Building
  @unpack x₁,x₂,x₃,x₄ = params
  dx = x₄[1]-x₃[1]
  dy = x₄[2]-x₃[2]
  building = gmsh.model.occ.addRectangle(x₃[1],x₃[2],0,dx,dy,tag=2)

  # Background Domain
  Lx = x₂[1]-x₁[1]
  Ly = x₂[2]-x₁[2]
  domain = gmsh.model.occ.addRectangle(x₁[1],x₁[2],0,Lx,Ly)

  # Subtract pieces from domain
  domain = gmsh.model.occ.cut((2,domain),building,-1,true,true)
  domain = domain[1][end][end]

  # Get entities
  fluid_domain = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]-0.1,-0.1,Lx+0.1,Ly+0.1,0.1,2)

  ### HERE!!!
  wall1 = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,H-0.1,-0.1,L+0.1,H+0.1,0.1,1)
  wall2 = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,-0.1,-0.1,L+0.1,0.1,0.1,1)
  inlet_point = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,-0.1,-0.1,0.1,H+0.1,0.1,0)
  inlet_line = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,-0.1,-0.1,0.1,H+0.1,0.1,1)
  outlet_point = gmsh.model.occ.getEntitiesInBoundingBox(L-0.1,-0.1,-0.1,L+0.1,H+0.1,0.1,0)
  outlet_line = gmsh.model.occ.getEntitiesInBoundingBox(L-0.1,-0.1,-0.1,L+0.1,H+0.1,0.1,1)
  monopile_point = gmsh.model.occ.getEntitiesInBoundingBox(Cx-R-0.1,Cy-R-0.1,-0.1,Cx+R+0.1,Cy+R+0.1,0.1,0)
  monopile_line = gmsh.model.occ.getEntitiesInBoundingBox(Cx-R-0.1,Cy-R-0.1,-0.1,Cx+R+0.1,Cy+R+0.1,0.1,1)


  # Get entity tags
  fluid_domain_tags = [ entity[2] for entity in fluid_domain]
  wall_tags = [ entity[2] for entity in wall1]
  append!(wall_tags , [ entity[2] for entity in wall2])
  inlet_point_tags = [ entity[2] for entity in inlet_point]
  inlet_line_tags = [ entity[2] for entity in inlet_line]
  outlet_point_tags = [ entity[2] for entity in outlet_point]
  outlet_line_tags = [ entity[2] for entity in outlet_line]
  monopile_point_tags = [ entity[2] for entity in monopile_point]
  monopile_line_tags = [ entity[2] for entity in monopile_line]

  println("fluid_domain", fluid_domain)

  # Physical group
  gmsh.model.addPhysicalGroup(0,inlet_point_tags,1)
  gmsh.model.addPhysicalGroup(1,inlet_line_tags,1)
  gmsh.model.addPhysicalGroup(0,outlet_point_tags,2)
  gmsh.model.addPhysicalGroup(1,outlet_line_tags,2)
  gmsh.model.addPhysicalGroup(1,wall_tags,3)
  gmsh.model.addPhysicalGroup(0,monopile_point_tags,4)
  gmsh.model.addPhysicalGroup(1,monopile_line_tags,4)
  pg1 = gmsh.model.addPhysicalGroup(2,fluid_domain_tags)#,5,"fluid")
  gmsh.model.setPhysicalName(2,pg1,"fluid")
  gmsh.model.setPhysicalName(0,1,"inlet")
  gmsh.model.setPhysicalName(1,1,"inlet")
  gmsh.model.setPhysicalName(0,2,"outlet")
  gmsh.model.setPhysicalName(1,2,"outlet")
  gmsh.model.setPhysicalName(1,3,"walls")
  gmsh.model.setPhysicalName(0,4,"monopile")
  gmsh.model.setPhysicalName(1,4,"monopile")


  # Synchronize
  gmsh.model.occ.synchronize()
  gmsh.model.geo.synchronize()

  # Define mesh size
  function meshSizeCallback(dim,tag,x,y,z,lc)
    if (Cx-R-dxLeft)<x<(Cx) && √((x-Cx)^2+(y-Cy)^2) < (R+dxLeft)
      dist = abs(√((x-Cx)^2+(y-Cy)^2) - R)/R
      return min(h_fine * (1 + decay_factor_left * (dist^decay_exponent_left)), h_coarse)
    elseif (Cx)<x<(Cx+R+dxRight) && (Cy-R-dyBottom)<y<(Cy+R+dyTop)
      dist = abs(√((x-Cx)^2+(y-Cy)^2) - R)/R
      return min(h_fine * (1 + decay_factor_right * (dist^decay_exponent_right)), h_coarse)
    else
      return h_coarse
    end
  end
  gmsh.model.mesh.setSizeCallback(meshSizeCallback)
  gmsh.model.mesh.generate()

  println(gmsh.model.getEntitiesForPhysicalGroup(2,5))

  # Finalize
  β2 = round(β;digits=2)
  α2 = round(α,digits=2)
  filename = "$L-$num_perforations-$β2-$α2.msh"
  meshes_path=ENV["PerforatedCylinder_MESHES"]
  mesh_file = joinpath(meshes_path,filename)
  gmsh.write(mesh_file)
  gmsh.finalize()

end
