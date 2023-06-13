using AdaptiveHierarchicalRegularBinning
using AdaptiveHierarchicalRegularBinning: make_tree
using AbstractTrees
using Test


@testset "Tree" begin
  dims = 2
  d, n = 8, 2^15
  V = rand(Float32, (dims==1 ? (n, d) : (d, n))...)
  R = sort(rand(UInt, n))
  I = collect(1:n)
  l = div(sizeof(eltype(R))*8, d)
  smlth=1
  scale = 1.0
  offset = fill(zero(eltype(V)), d)
  tree = make_tree(V, R, I, l, smlth, d, scale, offset; dims=dims)

  foreach(PreOrderDFS(tree)) do node
    lR = encpoints(node)
    lG = lR .>> (sizeof(eltype(lR))*8 - bitlen(node)*depth(node))
    @test all(lG .== nodevalue(node))
  end
end
