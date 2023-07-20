using AdaptiveHierarchicalRegularBinning
using AbstractTrees

# Setup
n = 100000
d = 20
dims = 2

dpt = 4
sml = 1000

X = randn(n, d)

if dims == 2
  X = X |> transpose |> collect
end

tree = regural_bin(UInt128, X, dpt, sml; dims=dims)

"""
mins, maxs = tightbounds(X; dims=1)
Find the tighest box that contains all points.

# Arguments
  - `X`: The cloud of points.

# Keyword Arguments
  - `dims`: dimension than enumerates the point coordinates. Defaults to 1.
"""
function tightbounds(X; dims=1)
  mins = minimum(X, dims=dims)
  maxs = maximum(X, dims=dims)
  return mins, maxs
end

# Callbacks
function leaf_cb(leaf)
  setcontext!(leaf, extrema(points(leaf); dims=leaddim(leaf))[:])
end

function node_cb(node)
  r = [ (Inf, -Inf) for _ in 1:bitlen(node) ]

  for child in children(node)
    ctx = getcontext(child)
    for (d, (m, M)) in enumerate(ctx)
      r[d] = (min(r[d][1], m), max(r[d][2], M))
    end
  end

  setcontext!(node, r)
end

# Evaluate
# Tree context -> Vector of size d, each element is a tuple of (min, max) per dim
applypostorder!(tree, leaf_cb, node_cb)

map(x -> x[1], getcontext(tree)) .- minimum(points(tree), dims=2)