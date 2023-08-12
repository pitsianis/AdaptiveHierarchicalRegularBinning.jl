using Test, AdaptiveHierarchicalRegularBinning, AbstractTrees


function test_dense_nodes(n, d, L, maxP)
  X = randn(d, n)
  tree = ahrb(X, L, maxP; QT=UInt128)

  @test all( i==n.nidx for (i,n) in enumerate(tree.info.nodes) if n.nidx != 0 )
  @test isperm(map(nindex, PostOrderDFS(tree)))
end


function test_configurations()
  maxP = 2^7
  @testset "Config ($n, $d, $L, $maxP)" for n in 10 .^ (4:6), d in 2:5, L in 10:5:25
    test_dense_nodes(n,d,L,maxP)
  end
end


@testset "Test Dense Nodes" begin
  test_configurations()
end
