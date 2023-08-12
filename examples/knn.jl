using NearestNeighbors, AdaptiveHierarchicalRegularBinning, AbstractTrees, Random
using Test

function knnsearch(C, Q, k)
  kdtree = KDTree(C; leafsize=10) #brutetree = BruteTree(C)
  return knn(kdtree, Q, k, true)
end

@inbounds @views function mergedstidx!(idx, dst, idx0, dst0, reidx)
  for i in 1:size(idx0, 1)
    j0, j1 = 1, 1
    while j0 <= length(idx0[1]) && j1 <= size(idx, 1)
      if dst0[i][j0] < dst[j1, i]
        # idx[j1+1:end, i] .= idx[j1:end-1, i]
        # dst[j1+1:end, i] .= dst[j1:end-1, i]
        for jj = size(idx, 1):-1:j1+1
          idx[jj, i] = idx[jj-1, i]
          dst[jj, i] = dst[jj-1, i]
        end
        idx[j1, i] = reidx[idx0[i][j0]]
        dst[j1, i] = dst0[i][j0]
        j0 += 1
        j1 += 1
      else
        j1 += 1
      end
    end
  end
end

mutable struct NNinfo
  maxdist::Float64
end
maxdist(node::SpatialTree) = getcontext(node).maxdist
maxdist(node::NNinfo) = node.maxdist

mutable struct GTinfo
  idx::Matrix{Int64}
  dst::Matrix{Float64}
end
GTinfo() = GTinfo(zeros(Int, 0, 0), zeros(Float64, 0, 0))
getidx(node::SpatialTree) = getglobalcontext(node).idx
getdst(node::SpatialTree) = getglobalcontext(node).dst

## 
@testset "elementary" begin

  # Random.seed!(0)
  n = 100_000
  k = 6

  for d = 1:4
    maxL = min(120 ÷ d, 25) 
    maxP = Int(ceil(sqrt(n)))

    X = randn(d, n)

    @inline prunepredicate(t, s) = qbox2boxdist(t, s) * t.info.scale > maxdist(t)

    @inline postconsolidate(t) = ~isleaf(t) ? getcontext(t).maxdist = maximum(maxdist.(t.info.context[cindices(t)])) : nothing

    @inline @views function processleafpair(t, s)
      C = points(s)
      Q = points(t)
      kk = min(k, size(C, 2))
      tpointidx = range(t)
      spointidx = range(s)
      idx0, dst0 = knnsearch(C, Q, kk)

      idx = getglobalcontext(tree).idx
      dst = getglobalcontext(tree).dst

      mergedstidx!(idx[:, tpointidx], dst[:, tpointidx], idx0, dst0, spointidx)
      getcontext(t).maxdist = maximum(dst[k, tpointidx])
    end

    println("AHRB (d = $d)")
    @time begin
      tree = AdaptiveHierarchicalRegularBinning.ahrb_fixed_length(X, maxL, maxP; QT=UInt128, ctxtype=NNinfo, gtctype=GTinfo)
      lookup = invperm(tree.info.perm)
    end
    println(tree)

    # run once to compile
    tree.info.context .= NNinfo.(Inf)
    getglobalcontext(tree).idx = zeros(Int, k, n)
    getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
    multilevelinteractions(tree, prunepredicate, processleafpair, postconsolidate)

    println("Brute force")
    idxs, dsts = @time knnsearch(points(tree), points(tree), k)
    idxs, dsts = hcat(idxs...), hcat(dsts...)

    # run again to time and check correctness
    println("DTT")
    getglobalcontext(tree).idx = zeros(Int, k, n)
    getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
    tree.info.context .= NNinfo.(Inf)
    @time multilevelinteractions(tree, prunepredicate, processleafpair, postconsolidate)

    idx = getglobalcontext(tree).idx
    dst = getglobalcontext(tree).dst

    @test all((X - points(tree)[:, lookup]) .== 0)
    @test all(idx .== idxs)
    @test dst ≈ dsts

    println("KD-Tree")
    kdtree = KDTree(X; leafsize=10)
    idxs, dsts = knn(kdtree, X, k, true)
    @time begin
      kdtree = KDTree(X; leafsize=10)
      idxs, dsts = knn(kdtree, X, k, true)
    end
  end
end