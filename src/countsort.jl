"""
$(TYPEDSIGNATURES)

Sequential countsort.
Inputs `Va`, `Ra`, `Ia` are not altered.

 - `Va`: The matrix to sort.
          Assuming `Va` is of size `(d, n)` then `d` is the number of dimensions and `n` is the number of elements.
 - `Ra`: The representation vector of `Va`. Must be of length `n`.
          Sorting will be executed based on this vector.
 - `Ia`: The current permutation of `Va`. Must be of length `n`.

 - `Vb`: A matrix to store the sorted version of `Va`. Must be of size `(d, n)`
 - `Rb`: A vector to store the sorted version of `Ra`. Must be of length `n`
 - `Ib`: A vector to store the resulting permutation. Must be of length `n`

  - `C`: A vector acting as a working space for the parallel countsort algorithm. Must be of length `maximum(Ra)`
"""
countsort_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractVector{<:Unsigned}} = @inbounds @views begin
  fill!(C, 0)

  for r in Ra
    C[r] += 1
  end

  cumsum!(C, C)

  for i in eachindex(Ra)
    j = C[Ra[i]]
    C[Ra[i]] -= 1

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

  - `C`: A matrix acting as a working space for the parallel countsort algorithm. Must be of size `(maximum(Ra), Threads.nthreads())`
"""
countsort_par_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractMatrix{<:Unsigned}} = @inbounds @views begin
  n  = length(Ra)
  np = Threads.nthreads()
  b  = cld(n, np)

  @inline local_range(k) = (k-1)*b + 1 : min(k*b, n)

  fill!(C, 0)

  Threads.@threads for k=1:np
    for i in local_range(k)
      C[Ra[i], k] += 1
    end
  end

  cumsum!(C, C, dims=2)
  cumsum!(C[:, end], C[:, end])
  C[2:end, 1:end-1] .+= C[1:end-1, end]

  Threads.@threads for k=1:np
    for i in local_range(k)
      j = C[Ra[i], k]
      C[Ra[i], k] -= 1

      Vb[:, j] .= Va[:, i]
      Rb[j]     = Ra[i]
      Ib[j]     = Ia[i]
    end
  end
end
