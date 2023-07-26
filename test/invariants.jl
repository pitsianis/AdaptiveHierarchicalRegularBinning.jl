using AdaptiveHierarchicalRegularBinning
using AdaptiveHierarchicalRegularBinning: make_tree
using AbstractTrees
using Test


@testset "datastructure invariants" begin
  n = 1_000
  maxP = 8
  for maxL = 6:2:10
    for d = 2:12

      X = randn(d, n)
      Xcopy = copy(X)
      tree = ahrb!(X, maxL, maxP; dims=2, QT=UInt128)

      # Original points are permuted
      @test Xcopy[:, tree.info.perm] == tree.info.points
      @test X === tree.info.points

      # all leaves have up to p points except the ones at the maxL level
      @test all(size(points(node), 2) <= maxP
                for node in PreOrderDFS(tree) if depth(node) < maxL && isleaf(node))

      # all leaves are leaves
      @test all(isleaf.(Leaves(tree)))

      # relationship of quantized and actual box centers and sides
      @test all(qbox(node) ≈ tree.info.scale * box(node) for node in PreOrderDFS(tree))

      # each node represents a contiquous group of points, groups are ordered in preorder DFS
      @test all(minimum(low.(children(node))) == low(node) &&
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

      @test all(isequal(getcontext(node), extrema(points(node); dims=leaddim(node))[:])
                for node in PreOrderDFS(tree))

      # root is only the first node ## BROKEN
      @test_broken isroot(tree)
    end
  end
end
