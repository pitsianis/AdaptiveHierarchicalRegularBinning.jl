# AdaptiveHierarchicalRegularBinning

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pitsianis.github.io/AdaptiveHierarchicalRegularBinning.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pitsianis.github.io/AdaptiveHierarchicalRegularBinning.jl/dev)
[![Build Status](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml?query=branch%3Amain)

Given a cloud of points in a potentially high dimensional metric space, we partition the points into a hierarchy of regular bins (hypercubes), by subdividing each hypercube adaptively up to a maximum L levels. A bin is not subdivided if it is at the maximum level `maxL` or contains up to `maxP` points.

*This package is under heavy development, and the API is not stable yet.*

## Example

Install the package with
```julia
] add https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl
```

and run the following code
```julia
using AdaptiveHierarchicalRegularBinning
n = 100_000
d = 20
X = rand(d, n)
Xcopy = copy(X)
maxL = 6
maxP = 32
tree = ahrb!(X, maxL, maxP; dims=2, QT=UInt128);

# Original points are permuted
  @assert Xcopy[:, tree.info.perm] == tree.info.points
  @assert X === tree.info.points

  # all leaves have up to p points except the ones at the maxL level
  @assert all(size(points(node), 2) <= maxP
            for node in PreOrderDFS(tree) if depth(node) < maxL && isleaf(node))

  # all leaves are leaves
  @assert all(isleaf.(Leaves(tree)))

  # relationship of quantized and actual box centers and sides
  @assert all(qbox(node) ≈ tree.info.scale * box(node) for node in PreOrderDFS(tree))

  # each node represents a contiquous group of points, groups are ordered in preorder DFS
  @assert all(minimum(low.(children(node))) == low(node) &&
            maximum(high.(children(node))) == high(node)
            for node in PreOrderDFS(tree) if !isleaf(node))

  # users can add application-specific information to the tree
  # axis-aligned bounding box (AABB)
  function boundingbox(node)
    if isleaf(node)
      setcontext!(node, extrema(points(node); dims=leaddim(node))[:])
    else
      bound(acc, ctx) = [(min(a[1], c[1]), max(a[2], c[2])) for (a, c) in zip(acc, ctx)]
      children(node) |>
      (C) -> mapreduce(getcontext, bound, C) |>
             (R) -> setcontext!(node, R)
    end
  end

  foreach(boundingbox, PostOrderDFS(tree))

  @assert all(isequal(getcontext(node), extrema(points(node); dims=leaddim(node))[:])
            for node in PreOrderDFS(tree))
```




