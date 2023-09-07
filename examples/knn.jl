using NearestNeighbors, AdaptiveHierarchicalRegularBinning, AbstractTrees, Random
using ThreadsX
using Test

function knnsearch(C, Q, k)
  kdtree = KDTree(C; leafsize=10) #brutetree = BruteTree(C)
  return knn(kdtree, Q, k, true)
end

@inbounds @views function mergedstidx!(
  idx,  # k-by-N nearest-neighbor indices, found so far
  dst,  # k-by-N nearest-neighbor distances, found so far
  idxtgt_new, # k-by-Ntgt nearest-neighbor indices, local to the source leaf, found in the current leaf2leaf interaction
  dsttgt_new, # k-by-Ntgt nearest-neighbor distances, found in the current leaf2leaf interaction
  reidx # M-by-1 index mapping from the local indices in the source leaf to the global index
)
  for i in 1:size(idxtgt_new, 1)
    j0, j1 = 1, 1
    while j0 <= length(idxtgt_new[1]) && j1 <= size(idx, 1)
      if dsttgt_new[i][j0] < dst[j1, i]
        for jj = size(idx, 1):-1:j1+1
          idx[jj, i] = idx[jj-1, i]
          dst[jj, i] = dst[jj-1, i]
        end
        idx[j1, i] = reidx[idxtgt_new[i][j0]]
        dst[j1, i] = dsttgt_new[i][j0]
        j0 += 1
        j1 += 1
      else
        j1 += 1
      end
    end
  end
end

function oracle_interaction_prc( 
  tree::SpatialTree{V,E,GTC,C}, 
  idx::Matrix{Ti}, 
  dst::Matrix{Td}
) where {V, E, GTC, C, Ti<:Integer, Td<:Real}

  n = size(idx, 2)

  leafidx = [range(i) for i in Leaves(tree)]
  # which leaves contain all neighbors of a point
  leavesperpoint = [mapreduce(x -> sum(x .>= first.(leafidx)), union, idx[:,i]) for i = 1:n]
  # which leaves contain all neighbors of the points of a leaf
  leavesperleaf = map(rng -> union(leavesperpoint[rng]...), leafidx)

  leaves = collect(Leaves(tree))
  L = length(leaves)
  # distance matrix between leaves
  Dmax = [AdaptiveHierarchicalRegularBinning.maxbox2boxdist(leaves[i], leaves[j]) for i in 1:L, j in 1:L]
  D = [box2boxdist(leaves[i], leaves[j]) for i in 1:L, j in 1:L]
  
  # what are the distances of the leaves that contain all neighbors of the points of a leaf
  distperleaf = [unique(Dmax[i, leavesperleaf[i]]) for i = 1:L]
  # how many leaves contain all neighbors
  oracle = sum(length.(leavesperleaf))
  
  # how many leaves have the same or smaller distance than the leaves that contain all neighbors
  reality_static = sum( 1:L ) do i
    sum(Dmax[i, :] .<= maximum(distperleaf[i]))
  end

  # how many leaves have the same or smaller distance than the furthest away kth neighbor of each target box
  reality_dynamic = sum( 1:L ) do i
    kth_largest_distance = maximum( dst[ end, range(leaves[i]) ] )
    sum( D[i, :] .<= kth_largest_distance )
  end

  return (reality_static, reality_dynamic, oracle)

end

mutable struct NNinfo
  maxdist::Float64
  num_node_interactions::Int64
  num_dist_calculations::Int64
end
NNinfo() = NNinfo(Inf, 0, 0)
maxdist(node::SpatialTree) = getcontext(node).maxdist
maxdist(node::NNinfo) = node.maxdist
num_node_interactions(node::SpatialTree) = getcontext(node).num_node_interactions
num_node_interactions(node::NNinfo) = node.num_node_interactions
num_dist_calculations(node::SpatialTree) = getcontext(node).num_dist_calculations
num_dist_calculations(node::NNinfo) = node.num_dist_calculations

mutable struct GTinfo
  idx::Matrix{Int64}
  dst::Matrix{Float64}
end
GTinfo() = GTinfo(zeros(Int, 0, 0), zeros(Float64, 0, 0))
getidx(node::SpatialTree) = getglobalcontext(node).idx
getdst(node::SpatialTree) = getglobalcontext(node).dst

## 

# Random.seed!(0)
n = 100_000
k = 6

trees_kept = []

for d = 4:4
  maxL = min(120 ÷ d, 25)
  maxP = Int(ceil(sqrt(n)))

  X = rand(d, n)

  @inline function prunepredicate(t, s)
    # only count box interactions that are not pruned
    cond = qbox2boxdist(t, s) * AdaptiveHierarchicalRegularBinning.scalar(t) > maxdist(t)
    getcontext(t).num_node_interactions += !cond    
    cond
  end

  @inline function postconsolidate(t)
    nidx = t.info.nodes[nindex(t)].pidx
    context = t.info.context
    while nidx != 0 # while not root
      oldvalue = context[nidx].maxdist
      newvalue = maximum(c -> context[c].maxdist, t.info.children[nidx])
      if newvalue < oldvalue
        context[nidx].maxdist = newvalue
        nidx = t.info.nodes[nidx].pidx
      else
        break
      end 
    end
  end

  @inline @views function processleafpair(t, s)
    C = points(s)
    Q = points(t)
    kk = min(k, size(C, 2))
    tpointidx = range(t)
    spointidx = range(s)

    # TODO: remove `size(points(t), 2)` factor; it can be multiplied at the
    # very end when we collect the dynamic distance-calculation statistics
    getcontext(t).num_dist_calculations += size(points(s), 2) * size(points(t), 2)

    # ------------------------------ 
    # TODO: use our own knnsearch here & on-the-fly merging
    idx0, dst0 = knnsearch(C, Q, kk)

    idx = getglobalcontext(tree).idx
    dst = getglobalcontext(tree).dst

    mergedstidx!(idx[:, tpointidx], dst[:, tpointidx], idx0, dst0, spointidx)
    # ------------------------------
    # TODO: mergedstidx should return maximum kth distance
    getcontext(t).maxdist = maximum(dst[k, tpointidx])
  end

  println("AHRB (d = $d)")
  @time begin
    tree = AdaptiveHierarchicalRegularBinning.ahrb_fixed_length(X, maxL, maxP; QT=UInt128, ctxtype=NNinfo, gtctype=GTinfo)
    lookup = invperm(tree.info.perm)
  end
  println(tree)

  # run once to compile
  tree.info.context .= NNinfo.()
  getglobalcontext(tree).idx = zeros(Int, k, n)
  getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
  multilevelinteractions(tree, tree, prunepredicate, processleafpair, postconsolidate)

  println("Brute force")
  idxs, dsts = @time knnsearch(points(tree), points(tree), k)
  idxs, dsts = hcat(idxs...), hcat(dsts...)

  # run again to time and check correctness
  println("DTT")
  getglobalcontext(tree).idx = zeros(Int, k, n)
  getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
  tree.info.context .= NNinfo.()
  @time multilevelinteractions(tree, tree, prunepredicate, processleafpair, postconsolidate)

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

  ## priority
  println("priority dtt")
  tree.info.context .= NNinfo.()
  getglobalcontext(tree).idx = zeros(Int, k, n)
  getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
  @time prioritymultilevelinteractions(tree, tree,
    box2boxdist, prunepredicate, processleafpair, postconsolidate)

  # read results
  idx = getglobalcontext(tree).idx
  dst = getglobalcontext(tree).dst
  # get golden results
  idxs, dsts = knnsearch(points(tree), points(tree), k)
  idxs, dsts = hcat(idxs...), hcat(dsts...)

  @test all(idx .== idxs)
  @test dst ≈ dsts

  ## parallel priority
  println("parallel priority dtt")
  tree.info.context .= NNinfo.()
  getglobalcontext(tree).idx = zeros(Int, k, n)
  getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
  @time ThreadsX.foreach(t -> prioritymultilevelinteractions(t, tree,
    box2boxdist, prunepredicate, processleafpair, postconsolidate), collect(Leaves(tree)))

  # read results
  idx = getglobalcontext(tree).idx
  dst = getglobalcontext(tree).dst
  # get golden results
  idxs, dsts = knnsearch(points(tree), points(tree), k)
  idxs, dsts = hcat(idxs...), hcat(dsts...)

  @test all(idx .== idxs)
  @test dst ≈ dsts

  push!(trees_kept, tree)

end
