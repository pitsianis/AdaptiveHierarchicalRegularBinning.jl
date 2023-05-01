using AdaptiveHierarchicalRegularBinning: make_tree
using Test


@testset "Tree" begin
  A = rand(UInt, 2^15)
  l = 8
  bitlen = 8
  @time tree = make_tree(A, l, bitlen)

  function iterate_nodes(f, tree, root=1)
    f(tree.info.nodes[root])

    for child in tree.info.children[root]
      iterate_nodes(f, tree, child)
    end
  end

  iterate_nodes(tree) do node
    G = @view A[node.lo:node.hi]
    G = G .>> (sizeof(eltype(G))*8 - bitlen*node.depth)
    @test all(G .== G[1])
  end
end
