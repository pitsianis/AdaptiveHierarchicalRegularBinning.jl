# I am olny saving this to keep it for further testing in the future
function point2boxDist(p, node)
  h = box(node)
  c = center(node)
  
  sqrt(sum(max.(0.0, max(c .- h .- p, p .- (c .+ h))).^2))

end
