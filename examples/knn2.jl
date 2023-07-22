using AdaptiveHierarchicalRegularBinning
using AbstractTrees

# Setup
n = 1_000_000
d = 10
dims = 2

dpt = 4
sml = 1000

X = randn(n, d)

if dims == 2
  X = X |> transpose |> collect
end

tree = regular_bin(UInt128, X, dpt, sml; dims=dims)

function tightbox(node)
  if isleaf(node)
    setcontext!(node, extrema(points(node); dims=leaddim(node))[:])
  else
    bound(acc, ctx) = [ (min(a[1], c[1]), max(a[2], c[2])) for (a, c) in zip(acc,ctx) ]
    children(node) |>
      (C) -> mapreduce(getcontext, bound, C) |>
      (R) -> setcontext!(node, R)
  end
end

foreach(tightbox, PostOrderDFS(tree))

map(x -> x[1], getcontext(tree)) .- minimum(points(tree), dims=dims)