using AdaptiveHierarchicalRegularBinning, AbstractTrees
using Test

@testset "DTT" begin
  n = 40_000
  for d = 1:5
    X = randn(d, n)
    tree = ahrb(X, Int(ceil(log(n))), Int(ceil(sqrt(n))); QT=UInt128, ctxtype=Vector{Tuple{Int64,Float64}})
    # println(tree)

    nl = leafcount(tree)
    lv = collect(Leaves(tree))
    D = [qbox2boxdist(l1, l2) for l1 in lv, l2 in lv]

    foreach(node -> setcontext!(node, Tuple{Int64,Float64}[]), PreOrderDFS(tree))
    foreach(node -> sizehint!(getcontext(node), nl), PreOrderDFS(tree))

    processpair!(t, s) = push!(getcontext(t), (nindex(s), qbox2boxdist(t, s)))

    multilevelinteractions(tree, tree, (t, s) -> false, processpair!, (t) -> nothing)
    @test all(map(x -> length(getcontext(x)), Leaves(tree)) .== nl)

    @test all(sort(D[:, 1]) .== sort([x[2] for x in getcontext(first(Leaves(tree)))]))
  end
end

@testset "DTT-Inf" begin
  n = 40_000
  for d = 1:5
    X = randn(d, n)
    tree = ahrb(X, Int(ceil(log(n))), Int(ceil(sqrt(n))); QT=UInt128, ctxtype=Vector{Tuple{Int64,Float64}})
    # println(tree)

    nl = leafcount(tree)
    lv = collect(Leaves(tree))
    D = [qbox2boxdistInf(l1, l2) for l1 in lv, l2 in lv]

    foreach(PreOrderDFS(tree)) do node
      setcontext!(node, Tuple{Int64,Float64}[])
      sizehint!(getcontext(node), nl)
    end

    processpair!(t, s) = push!(getcontext(t), (nindex(s), qbox2boxdistInf(t, s)))

    multilevelinteractions(tree, tree, (t, s) -> qbox2boxdistInf(t, s) > 1 / 8, processpair!, (t) -> nothing)

    @test all(filter(x -> x <= 1 / 8, sort(D[:, 1])) .== sort([x[2] for x in getcontext(first(Leaves(tree)))]))
  end
end