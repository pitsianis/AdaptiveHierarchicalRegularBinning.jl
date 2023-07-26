# the function neighbors(node1, node2) returns one of tree symbols
# :apart, :coincident, :adjacent
# 
# We use the Morton encodings and the level of the nodes to determine the outcome
# 
# Each dimension is checked separately.
# Two boxes are adjacent if they are not apart in any dimension.
# Two boxes are apart if they are apart in at least one dimension.   
# Otherwise they are coincident.
#
# The meaning of the symbols is as follows:
# :apart : the boxes are not adjacent in this dimension, at least one possibly 
#   empty box is in-between the two boxes 
# :coincident : the boxes are the same spatially in this dimension
# adjacent: the boxes are adjacent in this dimension, i.e. the end of one box
#   is the start of the other box


function neighborhood(a,b)

  d = abs.(qcenter(a) .- qcenter(b))
  r = qbox(a) + qbox(b)

  # println("r: $r d: $d")

  if any(d .> r)
    return :apart
  elseif all(d .== 0)
    return :coincident
  else
    return :adjacent
  end
end

function initusercontext!(t,v)
  setcontext!(t, v)
end

function adjacentSameLevel(t,s)
  hood = neighborhood(t,s)
  if hood == :apart
    return
  else
    if hood == :adjacent
      push!(getcontext(t), s.nidx)
    end
    # adjacent or coincident
    for t_child in children(t) 
      for s_child in children(s)
        adjacentSameLevel(t_child, s_child)
     end
    end
  end
end

X = randn(2,400)
tree = ahrb!(UInt128, X, 6, 25; dims=2)

foreach(node -> setcontext!(node,[]), PreOrderDFS(tree))

adjacentSameLevel(tree,tree)

for nd in Iterators.filter(isleaf, PreOrderDFS(tree)) # PreOrderDFS(tree)
  println("$(nd.nidx) : $([Int64(v) for v in getcontext(nd)])")
end
