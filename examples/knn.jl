using AdaptiveHierarchicalRegularBinning
using AbstractTrees


heap_down!(A, I, L) = @inbounds begin
  @inline get_max_idx(i) = let l=2i, r=2i+1, n=length(A)
    if l > n
      i
    elseif r > n || A[l] > A[r]
      l
    else
      r
    end
  end

  n = length(A)
  i = 1
  while i <= n
    max_idx = get_max_idx(i)
    A[max_idx] <= A[i] && break
    A[i], A[max_idx] = A[max_idx], A[i]
    I[i], I[max_idx] = I[max_idx], I[i]
    L[i], L[max_idx] = L[max_idx], L[i]
    i = max_idx
  end
end

heap_push!(A, I, L, b, j, d) = @inbounds begin
  b >= first(A) && return
  A[1], I[1], L[1] = b, j, d
  heap_down!(A, I, L)
end

function heap_append!(A, I, L, B, J, d)
  for (b, j) in zip(B, J)
    heap_push!(A, I, L, b, j, d)
  end
end


heap_sort!(A, I, L) = @inbounds begin
  n = length(A)
  rA = @view A[1:n]
  rI = @view I[1:n]
  rL = @view L[1:n]

  while true
    rA[1], rA[end] = rA[end], rA[1]
    rI[1], rI[end] = rI[end], rI[1]
    rL[1], rL[end] = rL[end], rL[1]

    n = length(rA)-1

    n <= 1 && break

    rA = @view A[1:n]
    rI = @view I[1:n]
    rL = @view L[1:n]

    heap_down!(rA, rI, rL)
  end
end

knn_impl!(tree, child_idx, point_idx, indices, distances, levels) = @inbounds begin
  nindex(tree) == 0 && return false

  dims = leaddim(tree)
  root = AbstractTrees.getroot(tree)

  corpus_idx = collect(range(tree))

  if child_idx != 0
    child = SpatialTree(TreeInfo(tree), child_idx)
    setdiff!(corpus_idx, range(child))
  end

  corpus = staticselectdim(points(root), Val(dims), corpus_idx)
  point  = staticselectdim(points(root), Val(dims), point_idx)

  dists = sqrt.(sum( (corpus .- point) .^ 2, dims=1)) # TODO: Other metrics

  ldistances = staticselectdim(distances, Val(dims), point_idx)
  lindices = staticselectdim(indices, Val(dims), point_idx)
  llevels  = staticselectdim(levels, Val(dims), point_idx)
  heap_append!(ldistances, lindices, llevels, dists, corpus_idx, depth(tree))


  ret = true
  if any(abs.(point .- center(tree)) .+ first(ldistances) .>= box(tree))
    ret = knn_impl!(AbstractTrees.parent(tree), nindex(tree), point_idx, indices, distances, levels)
  end

  if child_idx == 0
    heap_sort!(ldistances, lindices, llevels)
  end

  return ret
end


function knn(tree, query, k)
  points(tree) == query || throw("Corpus must be the same as query.")

  dims = leaddim(tree)
  n = size(query, leaddim(tree))
  kn = dims==1 ? (n,k) : (k,n)
  indices   = fill(-1,  kn)
  distances = fill(Inf, kn)
  levels    = fill(-1,  kn)

  @sync for leaf in Leaves(tree)
    Threads.@spawn begin
      for point_idx in range(leaf)
        knn_impl!(leaf, 0, point_idx, indices, distances, levels)
      end
    end
  end

  return indices, distances, levels
end