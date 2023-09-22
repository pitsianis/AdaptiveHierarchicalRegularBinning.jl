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