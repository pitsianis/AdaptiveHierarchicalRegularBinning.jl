using AdaptiveHierarchicalRegularBinning
using AdaptiveHierarchicalRegularBinning: make_tree
using AbstractTrees
using Test


function testinvariants(X, tree)

  maxL = tree.info.maxdepth
  maxP = tree.info.maxpoints
  # Original points are permuted
  @test X[:, tree.info.perm] == points(tree)

  @test isperm(tree.info.perm)

  # all leaves have up to p points except the ones at the maxL level
  @test all(size(points(node), 2) <= maxP
            for node in PreOrderDFS(tree) if depth(node) < maxL && isleaf(node))

  # all leaves are leaves
  @test all(isleaf.(Leaves(tree)))

  # relationship of quantized and actual box centers and sides
  @test all(qsidelength(node) ≈ sidelength(node) / tree.info.scale for node in PreOrderDFS(tree))

  # each node represents a contiquous group of points, groups are ordered in preorder DFS
  @test all(minimum(low.(children(node))) == low(node) &&
            maximum(high.(children(node))) == high(node)
            for node in PreOrderDFS(tree) if !isleaf(node))

  # root is only the first node
  @test nindex.(Iterators.filter(isroot, PreOrderDFS(tree))) == [1]

  # every point is enclosed within its box
  are_within_box = map( tree.info.nodes ) do ni
    node = SpatialTree(tree.info, nindex(ni))
    c = center(node); h = sidelength(node)
    lft = c .- h/2 .- eps(2.0^25); 
    rgt = c .+ h/2 .+ eps(2.0^25)
    minp = minimum( points(node); dims = 2 )
    maxp = maximum( points(node); dims = 2 )
    r = all( lft .<= minp .<= maxp .< rgt )
    r || println(hcat(lft, minp, maxp, rgt))
    r
  end

  @test all( are_within_box )

end

@testset "datastructure invariants [$method]" for method = ["fixed-length", "block-ecp"]
  n = 1_000
  maxP = 8
  for maxL = 6:2:10
    for d = 2:12

      X = randn() .* randn(d, n) .+ randn(d, 1)
      tree = ahrb(X, maxL, maxP; QT=UInt128, method)

      testinvariants(X, tree)
    end
  end
end

@testset "full tree datastructure invariants [$method]" for method = ["fixed-length", "block-ecp"]
  maxmaxL = 6
  for d = 1:3
    for maxL = 2:maxmaxL
      X = fulltree(maxL, d)
      tree = ahrb(X, maxmaxL, 2^0; QT=UInt32, method)

      testinvariants(X, tree)
    end
  end
end


@testset "tall tree datastructure invariants [$method]" for method = ["fixed-length", "block-ecp"]
  n = 1_000
  maxmaxL = 6
  for d = 1:3
    for maxL = 2:maxmaxL
      X = lineartreedata(n, d, maxL)
      tree = ahrb(X, maxmaxL, 2^0; QT=UInt32, method)

      testinvariants(X, tree)
    end
  end
end

