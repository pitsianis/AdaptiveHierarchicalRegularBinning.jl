# NOTE: After ping-pong sorting use P array to find the true result. `Ra[P .== 1] = Rb[P .== 1]`
"""
$(TYPEDSIGNATURES)

  - `Va`: The matrix to sort.
            Assuming `Va` is of size `(d, n)` then `d` is the number of dimensions and `n` is the number of elements.
  - `Ra`: The representation vector of `Va`. Must be of length `n`.
            Sorting will be executed based on this vector.
  - `Ia`: The current permutation of `Va`. Must be of length `n`.

  - `Vb`: A matrix to store the sorted version of `Va`. Must be of size `(d, n)`.
  - `Rb`: A vector to store the sorted version of `Ra`. Must be of length `n`.
  - `Ib`: A vector to store the resulting permutation. Must be of length `n`.

  - `P`: The ping-pong array.

  - `rsd`: The `RadixSortDetails` for this recursion.
"""
radixsort_seq_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, P::AbstractVector{UInt8}, rsd::RadixSortDetails) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}} = @inbounds @views begin

  is_deep(rsd)  && return
  is_small(rsd) && return

  # TODO: Type correctness between C, lo and hi
  C = alloc!(rsd, Int(bitmask(rsd)+1))

  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, CountSortDetails(rsd))
  P[low(rsd):high(rsd)] .= depth(rsd)%2

  for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : high(rsd)
    nrsd = next(rsd, nlo, nhi)

    radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, nrsd)
  end

  free!(rsd, C)
end


"""
$(TYPEDSIGNATURES)

  - `Va`: The matrix to sort.
            Assuming `Va` is of size `(d, n)` then `d` is the number of dimensions and `n` is the number of elements.
  - `Ra`: The representation vector of `Va`. Must be of length `n`.
            Sorting will be executed based on this vector.
  - `Ia`: The current permutation of `Va`. Must be of length `n`.

  - `Vb`: A matrix to store the sorted version of `Va`. Must be of size `(d, n)`.
  - `Rb`: A vector to store the sorted version of `Ra`. Must be of length `n`.
  - `Ib`: A vector to store the resulting permutation. Must be of length `n`.

  - `P`: The ping-pong array.

  - `rsd`: The `RadixSortDetails` for this recursion.
"""
radixsort_par_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, P::AbstractVector{UInt8}, rsd::RadixSortDetails) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}} = @inbounds @views begin

  C = alloc!(rsd, Int(bitmask(rsd)+1))
  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, CountSortDetails(rsd))
  P[low(rsd):high(rsd)] .= depth(rsd)%2

  @sync for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : high(rsd)

    nrsd = next(rsd, nlo, nhi)

    length(nrsd) <= 1 && continue

    if is_big(nrsd)
      Threads.@spawn radixsort_par_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nrsd)
    else
      Threads.@spawn radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nrsd)
    end
  end

  free!(rsd, C)
end

"""
$(TYPEDSIGNATURES)

  - `Va`: The matrix to sort.
            Assuming `Va` is of size `(d, n)` then `d` is the number of dimensions and `n` is the number of elements.
  - `Ra`: The representation vector of `Va`. Must be of length `n`.
            Sorting will be executed based on this vector.
  - `Ia`: The current permutation of `Va`. Must be of length `n`.

  - `Vb`: A matrix to store the sorted version of `Va`. Must be of size `(d, n)`.
  - `Rb`: A vector to store the sorted version of `Ra`. Must be of length `n`.
  - `Ib`: A vector to store the resulting permutation. Must be of length `n`.

  - `P`: The ping-pong array.

  - `rsd`: The `RadixSortDetails` for this recursion.
"""
radixsort_par_par_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, P::AbstractVector{UInt8}, rsd::RadixSortDetails) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}} = @inbounds @views begin

  Cm = alloc!(rsd, Int(bitmask(rsd)+1), Threads.nthreads())
  C  = alloc!(rsd, Int(bitmask(rsd)+1))
  countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, Cm, CountSortDetails(rsd))
  minimum!(reshape(C, :, 1), Cm)
  free!(rsd, Cm)

  P[low(rsd):high(rsd)] .= depth(rsd)%2

  @sync for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : high(rsd)
    nrsd = next(rsd, nlo, nhi)

    (is_small(nrsd) || is_deep(nrsd)) && continue

    if is_huge(nrsd)
      Threads.@spawn radixsort_par_par_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nrsd)
    elseif is_big(nrsd)
      Threads.@spawn radixsort_par_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nrsd)
    else
      Threads.@spawn radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nrsd)
    end
  end

  free!(rsd, C)
end


