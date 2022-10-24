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

  - `C`: A vector acting as a working space for the parallel countsort algorithm. Must be of length `maximum(Ra)+1`.

  - `lo`: Start index.
  - `hi`: End index.
"""
countsort_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC, lo::Unsigned, hi::Unsigned) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractVector{<:Unsigned}} = @inbounds @views begin
  fill!(C, 0)

  @inline to_index(x) = x + 1

  for i in lo:hi
    C[to_index(Ra[i])] += 1
  end

  cumsum!(C, C)
  C .+= lo - 1

  for i in lo:hi
    j = C[to_index(Ra[i])]
    C[to_index(Ra[i])] -= 1

    # TODO: `selectdim(Vb, dims, j) .= selectdim(Va, dims, i)`
    Vb[:, j] .= Va[:, i]
    Rb[j]     = Ra[i]
    Ib[j]     = Ia[i]
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

  - `Vb`: A matrix to store the sorted version of `Va`. Must be of size `(d, n)`
  - `Rb`: A vector to store the sorted version of `Ra`. Must be of length `n`
  - `Ib`: A vector to store the resulting permutation. Must be of length `n`

  - `C`: A matrix acting as a working space for the parallel countsort algorithm. Must be of size `(maximum(Ra)+1, Threads.nthreads())`

  - `lo`: Start index.
  - `hi`: End index.
"""
countsort_par_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC, lo::Unsigned, hi::Unsigned) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractMatrix{<:Unsigned}} = @inbounds @views begin
  n  = hi - lo + 1
  np = Threads.nthreads()
  b  = cld(n, np)

  @inline to_index(x) = x + 1
  @inline local_range(k) = (k-1)*b + lo : min(k*b, n) + lo - 1

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
      j = C[to_index(Ra[i]), k]
      C[to_index(Ra[i]), k] -= 1

      Vb[:, j] .= Va[:, i]
      Rb[j]     = Ra[i]
      Ib[j]     = Ia[i]
    end
  end
end
