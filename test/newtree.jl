using AdaptiveHierarchicalRegularBinning
using AbstractTrees
using Test

second((x, y)) = y

@testset "full tree" begin
  for maxL = 2:8
    r = 0.0:2.0^-(maxL - 1):(1.0-eps(1.0))
    P = [(x, y) for x = r, y = r]
    X = [first.(P[:])'; second.(P[:])']
    tree = ahrb!(X, 8, 2^0; dims=2, QT=UInt32)

    nnodes = treesize(tree)
    nleaves = length(collect(Leaves(tree)))

    @test nleaves == 2^(2 * (maxL - 1))
  end
end

