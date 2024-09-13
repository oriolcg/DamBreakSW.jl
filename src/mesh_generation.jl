"""
    struct Mesh_params

Parameters for the mesh generation. The following parameters are considered:
  - `filename::String`: Mesh file
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
  filename::String = "tmp_mesh.msh"      # Mesh file
  x₁::Vector = [-11.678,0.0]  # Lower left corner
  x₂::Vector = [32.0,1.4]     # Upper right corner
  x₃::Vector = [14.0,0.55]    # Lower left corner of the building (rectangular)
  x₄::Vector = [14.3,0.85]    # Upper right corner of the building (rectangular)
  H::Float64 = 0.5                       # Coarse element size
  h::Float64 = 0.1                       # Fine element size
  dxLeft::Float64 = 3.0                  # Refinement zone distance from object (left)
  dxRight::Float64 = 6.0                 # Refinement zone distance from object (right)
  dyTop::Float64 = 0.7                   # Refinement zone distance from object (top)
  dyBottom::Float64 = 0.7                # Refinement zone distance from object (bottom)
  decay_factor_left::Float64 = 2       # Decay factor for the left refinement zone
  decay_factor_right::Float64 = 2      # Decay factor for the right refinement zone
  decay_exponent_left::Float64 = 0.3     # Decay exponent for the left refinement zone
  decay_exponent_right::Float64 = 0.3    # Decay exponent for the right refinement zone
end

"""
    generate_mesh(params::Mesh_params)

Generate a mesh for the building problem.
"""
function generate_mesh(params::Mesh_params)

  # Initialize
  # gmsh.finalize()
  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  # gmsh.model.add("t1")

  # Building
  @unpack x₁,x₂,x₃,x₄ = params
  dx = x₄[1]-x₃[1]
  dy = x₄[2]-x₃[2]
  building = gmsh.model.occ.addRectangle(x₃[1],x₃[2],0,dx,dy)

  # Background Domain
  Lx = x₂[1]-x₁[1]
  Ly = x₂[2]-x₁[2]
  domain = gmsh.model.occ.addRectangle(x₁[1],x₁[2],0,Lx,Ly)

  # Subtract pieces from domain
  domain = gmsh.model.occ.cut((2,domain),(2,building),-1,true,true)
  println(domain)
  domain = domain[1][end][end]

  # Get entities
  fluid_domain = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]-0.1,-0.1,Lx+0.1,Ly+0.1,0.1,2)
  wall_points = gmsh.model.occ.getEntitiesInBoundingBox(x₃[1]-0.1,x₃[2]-0.1,-0.1,x₄[1]+0.1,x₄[2]+0.1,0.1,0)
  wall_lines = gmsh.model.occ.getEntitiesInBoundingBox(x₃[1]-0.1,x₃[2]-0.1,-0.1,x₄[1]+0.1,x₄[2]+0.1,0.1,1)
  inlet_points = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]-0.1,-0.1,x₁[1]+0.1,x₁[2]+Ly+0.1,0.1,0)
  inlet_line = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]-0.1,-0.1,x₁[1]+0.1,x₁[2]+Ly+0.1,0.1,1)
  outlet_points = gmsh.model.occ.getEntitiesInBoundingBox(x₂[1]-0.1,x₂[2]-Ly-0.1,-0.1,x₂[1]+0.1,x₂[2]+0.1,0.1,0)
  outlet_line = gmsh.model.occ.getEntitiesInBoundingBox(x₂[1]-0.1,x₂[2]-Ly-0.1,-0.1,x₂[1]+0.1,x₂[2]+0.1,0.1,1)
  side1_points = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]-0.1,-0.1,x₂[1]+0.1,x₂[2]-Ly+0.1,0.1,0)
  side1_line = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]-0.1,-0.1,x₂[1]+0.1,x₂[2]-Ly+0.1,0.1,1)
  side2_points = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]+Ly-0.1,-0.1,x₂[1]+0.1,x₂[2]+0.1,0.1,0)
  side2_line = gmsh.model.occ.getEntitiesInBoundingBox(x₁[1]-0.1,x₁[2]+Ly-0.1,-0.1,x₂[1]+0.1,x₂[2]+0.1,0.1,1)

  # Get entity tags
  fluid_domain_tags = [ entity[2] for entity in fluid_domain]
  wall_points_tags = [ entity[2] for entity in wall_points]
  wall_lines_tags = [ entity[2] for entity in wall_lines]
  inlet_points_tags = [ entity[2] for entity in inlet_points]
  inlet_line_tags = [ entity[2] for entity in inlet_line]
  outlet_points_tags = [ entity[2] for entity in outlet_points]
  outlet_line_tags = [ entity[2] for entity in outlet_line]
  side_points_tags = [ entity[2] for entity in side1_points]
  append!(side_points_tags,[ entity[2] for entity in side2_points])
  side_line_tags = [ entity[2] for entity in side1_line]
  append!(side_line_tags,[ entity[2] for entity in side2_line])
  println(fluid_domain_tags)
  println(wall_points_tags)
  println(wall_lines_tags)
  println(inlet_points_tags)
  println(inlet_line_tags)
  println(outlet_points_tags)
  println(outlet_line_tags)
  println(side_points_tags)
  println(side_line_tags)

  # Physical group
  gmsh.model.addPhysicalGroup(0,inlet_points_tags,1)
  gmsh.model.addPhysicalGroup(1,inlet_line_tags,1)
  gmsh.model.addPhysicalGroup(0,outlet_points_tags,2)
  gmsh.model.addPhysicalGroup(1,outlet_line_tags,2)
  gmsh.model.addPhysicalGroup(0,wall_points_tags,3)
  gmsh.model.addPhysicalGroup(1,wall_lines_tags,3)
  gmsh.model.addPhysicalGroup(0,side_points_tags,4)
  gmsh.model.addPhysicalGroup(1,side_line_tags,4)
  pg1 = gmsh.model.addPhysicalGroup(2,fluid_domain_tags)#,5,"fluid")
  gmsh.model.setPhysicalName(2,pg1,"fluid")
  gmsh.model.setPhysicalName(0,1,"inlet")
  gmsh.model.setPhysicalName(1,1,"inlet")
  gmsh.model.setPhysicalName(0,2,"outlet")
  gmsh.model.setPhysicalName(1,2,"outlet")
  gmsh.model.setPhysicalName(1,3,"walls")
  gmsh.model.setPhysicalName(0,4,"sides")
  gmsh.model.setPhysicalName(1,4,"sides")


  # Synchronize
  gmsh.model.occ.synchronize()
  gmsh.model.geo.synchronize()

  # Define mesh size
  @unpack H,h,dxLeft,dxRight,dyTop,dyBottom,decay_factor_left,decay_factor_right,decay_exponent_left,decay_exponent_right = params
  function meshSizeCallback(dim,tag,x,y,z,lc)
    Cx = (x₃[1]+x₄[1])/2
    Cy = (x₃[2]+x₄[2])/2
    R = (dx+dy)/2
    if dim==1 && tag in wall_lines_tags
      return h
    end
    if (Cx-dxLeft)<x<(Cx) && √((x-Cx)^2+(y-Cy)^2) < (dxLeft)
      dist = abs(√((x-Cx)^2+(y-Cy)^2))
      return min(h * (1 + decay_factor_left * (dist^decay_exponent_left)), H)
    elseif (Cx)<x<(Cx+dxRight) && (Cy-dyBottom)<y<(Cy+dyTop)
      dist = abs(√((x-Cx)^2+(y-Cy)^2))
      return min(h * (1 + decay_factor_right * (dist^decay_exponent_right)), H)
    else
      return H
    end
  end
  gmsh.model.mesh.setSizeCallback(meshSizeCallback)
  gmsh.model.mesh.generate()

  # println(gmsh.model.getEntitiesForPhysicalGroup(2,5))

  # Finalize
  @unpack filename = params
  mesh_file = datadir("meshes",filename)
  println("Mesh file: ",mesh_file)
  gmsh.write(mesh_file)
  gmsh.finalize()

end
