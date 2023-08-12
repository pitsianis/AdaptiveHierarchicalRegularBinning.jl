using DrWatson
@quickactivate "AHRB"

using AdaptiveHierarchicalRegularBinning, BenchmarkTools, Random, DataFrames, PythonPlot,
  ColorSchemes, Colors


"""
    rand_particles(num_particles::Int64)

Generate `num_particles` random particles using a distribution similar to a
disk-shaped galaxy.
"""
function rand_particles(num_particles::Int64)
  pos  = zeros(3, num_particles)
  vel  = zeros(3, num_particles)
  mass = zeros(num_particles)
  for i = 1:num_particles
    θ = 2π*rand()
    R = randexp()/ℯ
    z = (rand() > 0.5 ? 1 : -1) * (randexp()/(10ℯ))
    v = √(R*num_particles/1e6)
    p_mass = 1e4
    pos[:,i] .= [R*cos(θ), R*sin(θ), z]
    vel[:,i] .= [v*sin(θ), -v*cos(θ), 0.0]
    mass[i] = p_mass
  end
  return (pos, vel, mass)
end

"""
    show_particles(particles, θ = 30.0, ϕ = 30.0)

Using the `PythonPlot` package, display an array of particles in 3D with azimuthal
camera angle θ and elevation angle ϕ.
"""
function show_particles!(ax::Py, pos::Matrix{Float64}, θ::Float64 = 30.0, ϕ::Float64 = 30.0)
    
  ax.scatter(pos[1,:],pos[2,:],pos[3,:],s=0.4,c="black",alpha=0.4)
  ax.view_init(azim=θ, elev=ϕ)
  ax.set_xlim(-1.0, 1.0)
  ax.set_ylim(-1.0, 1.0)
  ax.set_zlim(-1.0, 1.0)

  ax.set_xticklabels("")
  ax.set_yticklabels("")
  ax.set_zticklabels("")

end

function show_particles!(axs::Matrix{Py}, pos::Matrix{Float64})
    
  map( x -> x.cla(), axs )
  map( x -> x.set_proj_type("persp", focal_length=0.2), axs )
  
  show_particles!(axs[1,1], pos, 30.0, 30.0)
  show_particles!(axs[1,2], pos, 10.0, 80.0)
  show_particles!(axs[2,1], pos, 80.0, 10.0)
  show_particles!(axs[2,2], pos, 60.0, 30.0)


end


getmass(node::SpatialTree) = getcontext(node)[:mass]
getcom(node::SpatialTree)  = getcontext(node)[:com]

"""
    grav_acc(mass, r; ϵ = 0.02)

Calculate the gravitation acceleration on a particle given the mass of another
object and the distance vector from particle to object. The usual gravitational
force is softened by a factor `ϵ` to prevent it from blowing up when the
distance is very small.
"""
function grav_acc(mass::Float64, r::Array{Float64,1}; ϵ::Float64 = 0.02)
    # Calculate gravitational acceleration, using node center of mass
    # and softening factor ϵ (to prevent blow up at d2 = 0)
    G::Float64 = 6.67430e-11
    ((G * mass) / ((sum(r.^2)+ϵ^2)^(3/2))) .* r
end

"""
    net_acc(particle, node, θ, acc_func)

Recursively calculate the net acceleration on a particle given the Barnes Hut
octree, a threshold `θ`, and the acceleration function `acc_func`.
"""
function net_acc(pos::Vector{Float64}, vel::Vector{Float64}, mass::Float64, node::SpatialTree, θ::Float64, all_masses::Vector{Float64})
  
  s = sidelength(node)     # Width of the cube that this node represents
  r = getcom(node) .- pos  # Vector from particle to node center of mass
  
  # Recusively calculates net acceleration on a particle given a particle tree
  if isleaf( node ) 
    # Base case 1: if leaf node, directly interact with all particles
    pp  = points(node)
    idx = range(node)

    acc = zeros( size(pos,1) )
    for i in axes(pp, 2)
      idx_orig = tree.info.perm[idx[i]]
      r  .= pp[:,i] .- pos
      acc .+= grav_acc( all_masses[idx_orig], rp )
    end

    return acc

  elseif s/√sum(r.^2) < θ  # expansion
    # Base case 2: if far-field, directly interact with box (expansion)
    return grav_acc( getmass(node), r )
  else
    # Recursive case, Barnes Hut method
    return sum( net_acc(pos, vel, mass, child, θ, all_masses) for child in children(node) )
  end

end

function simulation_step!(pos::Matrix{Float64}, vel::Matrix{Float64}, masses::Vector{Float64}, tree::SpatialTree, Δt::Float64, θ::Float64)
  # Calculate acceleration, and approximately advance velocity and position (using threading)
  for i in axes(pos,2)
    vel[:,i] += net_acc(pos[:,i], vel[:,i], masses[i], tree, θ, masses) * Δt
    pos[:,i] += vel[:,i] * Δt
  end
end

## run simulation

frames = Vector{Matrix{Float64}}()

pos,vel,masses = rand_particles(1000);
push!(frames, copy(pos))

for i in 1:1000
  tree = ahrb(pos, 10, 4; ctxtype = NamedTuple{(:com, :mass), Tuple{Vector{Float64}, Float64}} );

  foreach(PostOrderDFS(tree)) do node
    if isleaf(node)
      com  = points(node) * masses[ tree.info.perm[range(node)] ]
      mass = sum( @view( masses[ tree.info.perm[range(node)] ] ) )    
    else
      com = zeros(size(pos,1)); mass = 0.0;
      for child in children(node)
        com  .+= getcontext(child)[:com] .* getcontext(child)[:mass]
        mass  += getcontext(child)[:mass]
      end
    end
    com ./= mass
    setcontext!(node, (; com = com, mass = mass))
  end

  simulation_step!(pos, vel, masses, tree, 0.1, 0.5)
  push!(frames, copy(pos))
end

## animation

pyplot.rcParams["text.usetex"] = true

fig, axs = subplots(2, 2, layout="constrained", subplot_kw=pydict( Dict( "projection"=>"3d" ) ) )
fig.set_size_inches(8, 8)

pyanimation = pyimport("matplotlib.animation")

function animate(i)
  global axs
  global frames
  show_particles!( pyconvert(Matrix, axs), frames[i+1] )
  return
end

# myanim = pyanimation.FuncAnimation(fig, animate, frames=length(frames))
# myanim.save("bh-animation.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"], fps=20 )