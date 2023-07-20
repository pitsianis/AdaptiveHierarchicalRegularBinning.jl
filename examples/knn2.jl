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
