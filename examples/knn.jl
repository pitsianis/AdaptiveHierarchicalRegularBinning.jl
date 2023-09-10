using NearestNeighbors, AdaptiveHierarchicalRegularBinning, AbstractTrees, SparseArrays, Random
using ThreadsX
using Test

import AdaptiveHierarchicalRegularBinning: onthefly_block_knn!, heap_sort_inplace!

function knnsearch(C, Q, k)
  kdtree = KDTree(C; leafsize=10) #brutetree = BruteTree(C)
  return knn(kdtree, Q, k, true)
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

  @inline function processleafpair(t, s)
    C = points(s)
    Q = points(t)
    kk = min(k, size(C, 2))
    tpointidx = range(t)
    spointidx = range(s)

    # TODO: remove `size(points(t), 2)` factor; it can be multiplied at the
    # very end when we collect the dynamic distance-calculation statistics
    getcontext(t).num_dist_calculations += size(points(s), 2) * size(points(t), 2)
    
    idx = getglobalcontext(tree).idx
    dst = getglobalcontext(tree).dst

    idx_t = @view idx[:, tpointidx]
    dst_t = @view dst[:, tpointidx]
    
    # we should evaluate these once and store them (where to put these?)
    Q_nrmsq = vec(sum(abs2, Q, dims=1))
    C_nrmsq = vec(sum(abs2, C, dims=1))

    # must be sparse and we need to keep transpose of C (where to put these?)
    Ct = sparse( permutedims( C ) )
    Q  = sparse( Q )
    xb = zeros( size(C,2) )

    maxdist = onthefly_block_knn!( idx_t, dst_t, xb, Ct, Q, C_nrmsq, Q_nrmsq, spointidx )

    getcontext(t).maxdist = maxdist
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

  # final sorting
  @inbounds for i in axes(dst,2)
    heap_sort_inplace!( view(dst, :, i), view(idx, :, i) )
  end

  @test all((X - points(tree)[:, lookup]) .== 0)
  @test all(idx .== idxs)
  @test dst ≈ dsts.^2

  println("KD-Tree")
  kdtree = KDTree(X; leafsize=10)
  idxs, dsts = knn(kdtree, X, k, true)
  @time begin
    kdtree = KDTree(X; leafsize=10)
    idxs, dsts = knn(kdtree, X, k, true)
  end

  #= priority
  println("priority dtt")
  tree.info.context .= NNinfo.()
  getglobalcontext(tree).idx = zeros(Int, k, n)
  getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
  @time prioritymultilevelinteractions(tree, tree,
    box2boxdist, prunepredicate, processleafpair, postconsolidate)

  # read results
  idx = getglobalcontext(tree).idx
  dst = getglobalcontext(tree).dst

  # final sorting
  @inbounds for i in axes(dst,2)
    heap_sort_inplace!( view(dst, :, i), view(idx, :, i) )
  end

  # get golden results
  idxs, dsts = knnsearch(points(tree), points(tree), k)
  idxs, dsts = hcat(idxs...), hcat(dsts...)

  @test all(idx .== idxs)
  @test dst ≈ dsts.^2
  =#

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
  
  # final sorting
  @inbounds for i in axes(dst,2)
    heap_sort_inplace!( view(dst, :, i), view(idx, :, i) )
  end
  
  # get golden results
  idxs, dsts = knnsearch(points(tree), points(tree), k)
  idxs, dsts = hcat(idxs...), hcat(dsts...)

  @test all(idx .== idxs)
  @test dst ≈ dsts.^2

  ## specialized parallel priority
  println("specialized parallel priority dtt")
  tree.info.context .= NNinfo.()
  getglobalcontext(tree).idx = zeros(Int, k, n)
  getglobalcontext(tree).dst = ones(Float64, k, n) * Inf
  @time begin
    L = collect(Leaves(tree))
    ThreadsX.foreach(t -> processleafpair(t, t), L)
    ThreadsX.foreach(t -> specialprioritymultilevelinteractions(t, tree,
    box2boxdist, prunepredicate, processleafpair, postconsolidate), L)
  end

  # read results
  idx = getglobalcontext(tree).idx
  dst = getglobalcontext(tree).dst
  
  # final sorting
  @inbounds for i in axes(dst,2)
    heap_sort_inplace!( view(dst, :, i), view(idx, :, i) )
  end
  
  # get golden results
  idxs, dsts = knnsearch(points(tree), points(tree), k)
  idxs, dsts = hcat(idxs...), hcat(dsts...)

  @test all(idx .== idxs)
  @test dst ≈ dsts.^2

  push!(trees_kept, tree)

end
