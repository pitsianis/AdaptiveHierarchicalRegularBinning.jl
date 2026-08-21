using AdaptiveHierarchicalRegularBinning
using AbstractTrees
using Test


@testset "full tree" begin
  maxmaxL = 6
  for d = 1:3
    for maxL = 2:maxmaxL
      X = fulltree(maxL, d)
      tree = ahrb(X, maxmaxL, 2^0; QT=UInt32)

      nnodes = treesize(tree)
      nleaves = length(collect(Leaves(tree)))

      @test treeheight(tree) == maxL
      @test nleaves == 2^(d * maxL)
      @test nnodes == sum((2^d) .^ l for l = 0:maxL)
    end
  end
end

@testset "tall tree" begin
  n = 1_000
  maxmaxL = 6
  for d = 1:3
    for maxL = 2:maxmaxL
      X = lineartreedata(n, d, maxL)
      tree = ahrb(X, maxmaxL, 2^0; QT=UInt32)

      # println("d = $d, maxL = $maxL breadth = $(treebreadth(tree))")

      @test treeheight(tree) >= maxL
    end
  end
end

@testset "node context" begin
  tree = AdaptiveHierarchicalRegularBinning.ahrb(rand(8,10000), 12, 128; ctxtype=Int64)
  tree.info.context .= 1:treesize(tree)
  @test all(i -> getcontext(tree,i) == i, 1:treesize(tree) )
end

@testset "original indices" begin
  X = rand(3, 100)
  tree = ahrb(X, 5, 4; QT=UInt32)

  @test X[:, original_indices(tree)] == points(tree)
  @test isperm(original_indices(tree))

  node = first(Leaves(tree))
  indices = original_indices(node)
  @test X[:, indices] == points(node)
  indices[1] = 0
  @test X[:, original_indices(node)] == points(node)
end

@testset "leaf ranges" begin
  X = rand(3, 100)
  tree = ahrb(X, 5, 4; QT=UInt32)
  leaves = collect(Leaves(tree))

  @test collect(leafranges(tree)) == range.(leaves)
  original = collect(leafranges(tree; original=true))
  @test all(X[:, indices] == points(node) for (node, indices) in zip(leaves, original))
end

@testset "occupancy counts" begin
  X = fulltree(3, 2)
  tree = ahrb(X, 3, 1; QT=UInt32)

  @test occupancy_counts(tree) == [4, 16, 64, 0]
  @test occupancy_counts(tree; levels=[0, 2, 5]) == [1, 16, 0]
end