using SparseArrays

import AdaptiveHierarchicalRegularBinning: onthefly_block_knn!, heap_sort_inplace!

function evaluate(x::AbstractVector, y::AbstractVector)
  return sum(abs2.(x .- y))
end

function evaluate(x::AbstractVector, y::AbstractMatrix)
  return vec( sum(abs2.(x .- y); dims = 1) )
end

function ground_truth_knn(C::AbstractMatrix, Q::AbstractMatrix, k::Integer)
  n = size(Q, 2)
  idxs = zeros(Int, k, n)
  dists = zeros(eltype(Q), k, n)
  for i in 1:n
    d = evaluate(Q[:,i], C)
    ii = sortperm(d)[1:k]
    dists[:,i], idxs[:,i] = d[ii], ii
  end
  return idxs, dists
end

function onthefly_knn(C::AbstractMatrix, Q::AbstractMatrix, 
  C_nrmsq::AbstractVector, Q_nrmsq::AbstractVector, k::Integer)

  idxs  = zeros(Int, k, size(Q,2))
  dists = fill(typemax(eltype(Q)), k, size(Q,2))

  # partition the corpus set into equal-sized blocks
  @inbounds for blk_rng::UnitRange{Int64} in Iterators.partition( 1:size(C,2), 1000 )
      
    Cbt = permutedims( C[ :, blk_rng ] )
    Cb_nrmsq = @view C_nrmsq[ blk_rng ]
    xb = zeros( size(Cbt,1) )

    maxdist = onthefly_block_knn!( idxs, dists, xb, Cbt, Q, Cb_nrmsq, Q_nrmsq, blk_rng )

  end

  @inbounds for i in axes(dists,2)
    heap_sort_inplace!( view(dists, :, i), view(idxs, :, i) )
  end

  return idxs, dists

end

@testset "On-the-fly kNN" begin

  d  = 3
  m  = 10000
  k  = 4
  dt = Float64

  X = sparse( rand(dt, d, m) )

  X_nrmsq = vec(sum(abs2, X , dims=1))

  Q = X[ :, 100:200 ]
  Q_nrmsq = X_nrmsq[ 100:200 ]

  idxs, dists = onthefly_knn( X, Q, X_nrmsq, Q_nrmsq, k )

  idxs_gt, dists_gt = ground_truth_knn( X, Q, k )
  
  @test idxs == idxs_gt
  @test dists ≈ dists_gt

end
