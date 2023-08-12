"""
box2boxdist(c1, h1, c2, h2)

box2boxdist(node1, node2)

Compute the distance between two boxes.
c denotes the center of the box, h denotes the half-width of the box.

"""
function box2boxdist(c1, h1, c2, h2)
  sqrt(sum(max.(c2 .- h2 .- (c1 .+ h1), c1 .- h1 .- (c2 .+ h2), 0.0).^2))
end
function box2boxdist(node1, node2)
  c1 = center(node1)
  c2 = center(node2) 
  h = (sidelength(node1) + sidelength(node2))/2
  sqrt(sum(max.(c2 .- c1 .- h, c1 .- c2 .- h, 0.0).^2))
end


"""
point2boxdist(p, c, h)

point2boxdist(p, node)

Compute the distance between a point and a box.
p denotes the point.
c denotes the center of the box, h denotes the half-width of the box.
"""
function point2boxdist(p, c, h)
  sqrt(sum(max.(c .- h .- p, p .- (c .+ h), 0.0).^2))
end
function point2boxdist(p, node)
  c, h = center(node), sidelength(node)/2
  sqrt(sum(max.(c .- h .- p, p .- (c .+ h), 0.0).^2))
end

function qbox2boxdist(node1, node2)
  c1 = qcenter(node1)
  c2 = qcenter(node2) 
  h = (qsidelength(node1) + qsidelength(node2))/2
  # sqrt(sum(max.(c2 .- c1 .- h, c1 .- c2 .- h, 0.0).^2))
  # the following allocates less memory
  sm = zero( eltype(c1) )
  
  @inbounds for i in eachindex(c1)
    mx = zero(eltype(c1))
    mx = max(mx, c2[i] - c1[i] - h)
    mx = max(mx, c1[i] - c2[i] - h)
    sm += mx^2
  end

  return sqrt(sm)
end

function qbox2boxdistInf(node1, node2)
  c1 = qcenter(node1)
  c2 = qcenter(node2) 
  h = (qsidelength(node1) + qsidelength(node2))/2
  # sqrt(sum(max.(c2 .- c1 .- h, c1 .- c2 .- h, 0.0).^2))
  # the following allocates less memory
  sm = zero( eltype(c1) )
  
  @inbounds for i in eachindex(c1)
    mx = zero(eltype(c1))
    mx = max(mx, c2[i] - c1[i] - h)
    mx = max(mx, c1[i] - c2[i] - h)
    sm = max(sm,mx)
  end

  return sm
end
