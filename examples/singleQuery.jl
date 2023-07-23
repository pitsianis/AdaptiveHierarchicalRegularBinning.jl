## A snipet to show how to perform a single query
using AdaptiveHierarchicalRegularBinning, AbstractTrees

±(dx) = (x) -> (x+dx, x-dx)
square(x) = x .^ 2
dist(x, y) = sqrt(sum(square ∘ splat(-), zip(x, y)))

function point2boxDist(p, node)
  h = box(node)
  c = center(node)
  if all(abs.(p .- c) .<= h)
    return 0.0
  else
    return sqrt(sum(splat(min) ∘ square ∘ ±(h) ∘ splat(-), zip(c, p)))
  end
end

function predicate(node, query, r)
  return r >= point2boxDist(query, node)
end

function searchTree(tree, query, r=Inf)
  println("Visiting node ", tree.nidx, " with half-side ", box(tree))
  if isleaf(tree)
    dists = map((p) -> dist(p, query), eachslice(points(tree); dims=leaddim(tree)))
    k = argmin(dists)
    return dists[k], range(tree)[k]
  else
    cc = collect(children(tree))
    dd = [point2boxDist(query, node) for node in cc]
    visitorder = sortperm(dd)
    j = -1
    for i in visitorder
      if predicate(cc[i], query, r)
        r, j = searchTree(cc[i], query, r)
      end
    end
    return r, j
  end
end

d = 5; n = 40000
X = rand(d, n)
tree = regular_bin(UInt128, X, 6, 2^5; dims=2)
q = rand(d)
r = Inf
r, j = searchTree(tree, q, r)
i = tree.info.perm[j]

minarg((sqrt.(sum((points(tree) .- q).^2,dims=1)))...)

fig, axs = subplots(layout="constrained", figsize=(10,10))
ax.cla(); plottree(ax, tree)
ax.scatter( X[1,:], X[2,:], color="black", s=0.1)
display( fig )