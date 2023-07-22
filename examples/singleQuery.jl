## A snipet to show how to perform a single query
using AdaptiveHierarchicalRegularBinning, AbstractTrees

mutable struct Query
  coords::Vector{Float64}
  rho::Float64
end

function point2boxDist(p, node)
  h = box(node)
  c = center(node)
  if all(abs.(p .- c) .<= h)
    return 0.0
  else
    # compute the distance to the closest corner of the box (c +/- h)
    d = 0.0
    for i = 1:length(p)
      d += min(c[i]-h-p[i], c[i]+h-p[i])^2
    end
    d = sqrt(d)
  end
  return d
end

function predicate(node, query::Query)
  return query.rho >= point2boxDist(query.coords, node)
end

function updatesearch!(query, points)
  for p in points
    d = sqrt(sum((p .- query.coords).^2))
    if d < query.rho
      query.rho = d
      println("New closest point with distance ", d)
    end
  end
end

function searchTree(tree, query)
  println("Visiting node ", tree.nidx, " with half-side ", box(tree))
  if isleaf(tree)
    updatesearch!(query, points(tree))
  else
    cc = collect(children(tree))
    dd = [point2boxDist(query.coords, node) for node in cc]
    visitorder = sortperm(dd)
    for i in visitorder
      if predicate(cc[i], query)
        searchTree(cc[i], query)
      end
    end
  end
end

d = 2; n = 400
X = rand(d, n)
tree = regular_bin(UInt128, X, 6, 2^5; dims=2)
q = Query(rand(d),Inf)
searchTree(tree, q)

minimum(sqrt.(sum((points(tree) .- q.coords).^2,dims=1)))