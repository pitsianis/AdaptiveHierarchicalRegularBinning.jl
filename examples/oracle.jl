using NearestNeighbors, AdaptiveHierarchicalRegularBinning, AbstractTrees, Random

function knnsearch(C, Q, k)
  kdtree = brutetree = BruteTree(C) # KDTree(C; leafsize=10) # 
  return knn(kdtree, Q, k, true)
end

# Random.seed!(0)
n = 30_000
k = 6

for d = 2:20

  X = randn(d, n)
  tree = ahrb(X)

  # do knn with the permuted points to avoid renamings
  idx, dst = knnsearch(points(tree), points(tree), k)

  leafidx = [range(i) for i in Leaves(tree)]
  # which leaves contain all neighbors of a point
  leavesperpoint = [mapreduce(x -> sum(x .>= first.(leafidx)), union, idx[i]) for i = 1:n]
  # which leaves contain all neighbors of the points of a leaf
  leavesperleaf = map(rng -> union(leavesperpoint[rng]...), leafidx)

  leaves = collect(Leaves(tree))
  L = length(leaves)
  # distance matrix between leaves
  D = [box2boxdist(leaves[i], leaves[j]) for i in 1:L, j in 1:L]
  # what are the distances of the leaves that contain all neighbors of the points of a leaf
  distperleaf = [unique(D[i, leavesperleaf[i]]) for i = 1:L]
  # how many leaves contain all neighbors
  oracle = sum(length.(leavesperleaf))
  # how many leaves have the same or smaller distance than the leaves that contain all neighbors
  reality = sum([sum(D[i, :] .<= maximum(distperleaf[i])) for i = 1:L])

  println("$d : $((oracle, reality) ./ L^2)")
end