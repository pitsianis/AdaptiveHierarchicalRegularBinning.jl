# TODO: comments
"""
$(TYPEDSIGNATURES)

Sequential countsort.
Inputs `Va`, `Ra`, `Ia` are not altered.

  - `Va`: The matrix to sort.
            Assuming `Va` is of size `(d, n)` then `d` is the number of dimensions and `n` is the number of elements.
  - `Ra`: The representation vector of `Va`. Must be of length `n`.
            Sorting will be executed based on this vector.
  - `Ia`: The current permutation of `Va`. Must be of length `n`.

  - `Vb`: A matrix to store the sorted version of `Va`. Must be of size `(d, n)`.
  - `Rb`: A vector to store the sorted version of `Ra`. Must be of length `n`.
  - `Ib`: A vector to store the resulting permutation. Must be of length `n`.

  - `C`: A vector acting as a working space for the parallel countsort algorithm. Must be of length `256`.

  - `lo`: Start index.
  - `hi`: End index.

  - `nbyte`: The n-th most significant byte.
"""
countsort_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC, csd::CountSortDetails) where {TV<:AbstractArray, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractVector{<:Unsigned}} = @inbounds @views begin
  lo = low(csd)
  hi = high(csd)

  to_index(x) = radixsel(csd, x) + 1

  fill!(C, 0)
  for i in lo:hi
    C[to_index(Ra[i])] += 1
  end

  cumsum!(C, C)
  C .+= lo - 1

  for i in lo:hi
    ci = to_index(Ra[i])

    j = C[ci]
    C[ci] -= 1

    # selectdim(Vb, csd, j) .= selectdim(Va, csd, i)
    Rb[j] = Ra[i]
    # Ib[j] = Ia[i]
  end
end


"""
$(TYPEDSIGNATURES)

Multithreaded countsort.
Inputs `Va`, `Ra`, `Ia` are not altered.

  - `Va`: The matrix to sort.
            Assuming `Va` is of size `(d, n)` then `d` is the number of dimensions and `n` is the number of elements.
  - `Ra`: The representation vector of `Va`. Must be of length `n`.
            Sorting will be executed based on this vector.
  - `Ia`: The current permutation of `Va`. Must be of length `n`.

  - `Vb`: A matrix to store the sorted version of `Va`. Must be of size `(d, n)`.
  - `Rb`: A vector to store the sorted version of `Ra`. Must be of length `n`.
  - `Ib`: A vector to store the resulting permutation. Must be of length `n`.

  - `C`: A matrix acting as a working space for the parallel countsort algorithm. Must be of size `(256, Threads.nthreads())`.

  - `lo`: Start index.
  - `hi`: End index.

  - `f`: The transformation to apply.
"""
countsort_par_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC, csd::CountSortDetails) where {TV<:AbstractArray, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractMatrix{<:Unsigned}} = @inbounds @views begin
  lo = low(csd)
  hi = high(csd)

  n  = hi - lo + 1
  np = Threads.nthreads()
  b  = cld(n, np)

  to_index(x) = radixsel(csd, x) + 1
  local_range(k) = (k-1)*b + lo : min(k*b, n) + lo - 1

  fill!(C, 0)
  Threads.@threads for k=1:np
    for i in local_range(k)
      C[to_index(Ra[i]), k] += 1
    end
  end

  cumsum!(C, C, dims=2)
  cumsum!(C[:, end], C[:, end])
  C[2:end, 1:end-1] .+= C[1:end-1, end]
  C .+= lo - 1

  Threads.@threads for k=1:np
    for i in local_range(k)
      ci = to_index(Ra[i])

      j = C[ci, k]
      C[ci, k] -= 1

      # selectdim(Vb, csd, j) .= selectdim(Va, csd, i)
      Rb[j] = Ra[i]
      # Ib[j] = Ia[i]
    end
  end
end
