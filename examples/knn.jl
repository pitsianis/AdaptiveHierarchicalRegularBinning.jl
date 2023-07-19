using AdaptiveHierarchicalRegularBinning
using AbstractTrees


heap_down!(A, I, i=1) = @inbounds begin
  @inline get_max_idx(i) = @inbounds let l=2i, r=2i+1, n=length(A)
    if l > n
      i
    elseif r > n || A[l] > A[r]
      l
    else
      r
    end
  end

  n = length(A)
  while i <= n
    max_idx = get_max_idx(i)
    A[max_idx] <= A[i] && break
    A[i], A[max_idx] = A[max_idx], A[i]
    I[i], I[max_idx] = I[max_idx], I[i]
    i = max_idx
  end
end

heap_push!(A, I, b, j) = @inbounds begin
  b >= first(A) && return
  A[1], I[1] = b, j
  heap_down!(A, I)
end

function heap_append!(A, I, B, J, r)
  for (b, j) in zip(B, J)
    b >= r*r && continue
    heap_push!(A, I, b, j)
  end
end


heap_sort!(A, I) = @inbounds begin
  n = length(A)
  rA = @view A[1:n]
  rI = @view I[1:n]

  while true
    rA[1], rA[end] = rA[end], rA[1]
    rI[1], rI[end] = rI[end], rI[1]

    n = length(rA)-1

    n <= 1 && break

    rA = @view A[1:n]
    rI = @view I[1:n]

    heap_down!(rA, rI)
  end
end

heapify!(A, I) = @inbounds begin
  # Start from the last non-leaf node and move up to the root
  for i in (div(length(A), 2)):-1:1
    heap_down!(A, I, i)
  end
end


merge_heaps!(A, I, B, J) = begin
  C = Array{eltype(A)}(undef, length(A)+length(B))
  K = Array{eltype(I)}(undef, length(A)+length(B))
  i = 1  # Pointer for the existing heap
  j = 1  # Pointer for the heapified unsorted array

  # Merge the two heaps until one of them is exhausted
  while i <= length(A) && j <= length(B)
    if A[i] > B[j]
      C[i+j-1] = A[i]
      K[i+j-1] = I[i]
      i += 1
    else
      C[i+j-1] = B[j]
      K[i+j-1] = J[j]
      j += 1
    end
  end

  while i <= length(A)
    C[i+j-1] = A[i]
    K[i+j-1] = I[i]
    i += 1
  end

  while j <= length(B)
    C[i+j-1] = B[j]
    K[i+j-1] = J[j]
    j += 1
  end


  # TODO: Extract sub heap
  # A .= C[end-length(A)+1:end]
  # I .= K[end-length(I)+1:end]
  # heapify!(A, I)
end

knn_impl!(tree, child_idx, point_idx, indices, distances, r) = @inbounds begin
  nindex(tree) == 0 && return
  length(point_idx) == 0 && return

  dims = leaddim(tree)
  root = AbstractTrees.getroot(tree)

  corpus_idx = collect(range(tree))

  if child_idx != 0
    child = SpatialTree(TreeInfo(tree), child_idx)
    setdiff!(corpus_idx, range(child))
  end

  corpus = staticselectdim(points(root), Val(dims), corpus_idx)
  point  = staticselectdim(points(root), Val(dims), point_idx)

  c2 = sum((x) -> x^2, corpus; dims=dims==1 ? 2 : 1)
  p2 = sum((x) -> x^2, point;  dims=dims==1 ? 2 : 1)
  dists = abs.(transpose(c2) .- 2*transpose(corpus)*point .+ p2)

  i = 1 # point_idx index
  c = 1 # global counter
  while i <= length(point_idx)
    pnt_idx = point_idx[i]
    lpoint = staticselectdim(point, Val(dims), c)
    ldistances = staticselectdim(distances, Val(dims), pnt_idx)
    lindices = staticselectdim(indices, Val(dims), pnt_idx)
    ldists = staticselectdim(dists, Val(dims), c)

    heap_append!(ldistances, lindices, ldists, corpus_idx, r)

    # TODO: optimize heap_append!, profiler: 77% on all nn. O(mlog(n))
    # heapify!(ldists, corpus_idx) O(m)
    # merge_heaps!(ldistances, lindices, ldists, corpus_idx) O(m+n)

    effective_r = min(r, sqrt(first(ldistances)))

    if all(abs.(lpoint .- AdaptiveHierarchicalRegularBinning.center(tree)) .+ effective_r .<= box(tree))
      deleteat!(point_idx, i)
    else
      i += 1
    end
    c += 1
  end

  return knn_impl!(AbstractTrees.parent(tree), nindex(tree), point_idx, indices, distances, r)
end


function knn(tree, query; k=length(tree), r=Inf, sortres=false)
  points(tree) == query || throw("Corpus must be the same as query.")

  dims = leaddim(tree)
  n = size(query, leaddim(tree))
  kn = dims==1 ? (n,k) : (k,n)
  indices   = fill(-1,  kn)
  distances = fill(Inf, kn)

  @sync for leaf in Leaves(tree)
    Threads.@spawn knn_impl!(leaf, 0, collect(range(leaf)), indices, distances, r)
  end

  distances .= sqrt.(distances)

  if sortres
    nt = Threads.nthreads()
    kb = cld(n, nt)
    Threads.@threads for k in 1:nt
      for i in ((k-1)*kb + 1) : min(n, kb*k)
        ldistances = staticselectdim(distances, Val(dims), i)
        lindices = staticselectdim(indices, Val(dims), i)
        heap_sort!(ldistances, lindices)
      end
    end
  end

  return indices, distances
end