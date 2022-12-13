# NOTE: After ping-pong sorting use P array to find the true result. `Ra[P] .= Rb[P]`
"""
$(SIGNATURES)

Sequential-to-Sequential radixsort.

# Arguments
  - `Va`: The cloud of points.
  - `Ra`: The morton transformation for each point of `Va`
  - `Ia`: The current permutation of `Va`.

  - `Vb`: An auxiliary matrix to store the partially sorted version of `Va`.
  - `Rb`: An auxiliary vector to store the partially sorted version of `Ra`.
  - `Ib`: An auxiliary vector to store the resulting permutation.

  - `P`: A ping-pong matrix that denotes where the latest results are.

  - `rsd`: The `RadixSortDetails` object used to configure the algorithm.
"""
radixsort_seq_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, P::AbstractVector{Bool}, rsd::RadixSortDetails) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}} = @inbounds @views begin
  lo = low(rsd)
  hi = high(rsd)

  C = alloc!(rsd, Int(bitmask(rsd))+1)

  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, CountSortDetails(rsd))
  P[lo:hi] .= !P[lo]

  nC = length(C)
  for j in 1:nC
    nlo = C[j]+1
    nhi = j != nC ? C[j+1] : hi
    nrsd = next(rsd, nlo, nhi)

    (is_small(nrsd) || is_deep(nrsd)) && continue

    radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, nrsd)
  end

  free!(rsd, C)
end

#TODO: FIX bug with Threads.@spawn
"""
$(SIGNATURES)

Parallel-to-Sequential radixsort.

# Arguments
  - `Va`: The cloud of points.
  - `Ra`: The morton transformation for each point of `Va`
  - `Ia`: The current permutation of `Va`.

  - `Vb`: An auxiliary matrix to store the partially sorted version of `Va`.
  - `Rb`: An auxiliary vector to store the partially sorted version of `Ra`.
  - `Ib`: An auxiliary vector to store the resulting permutation.

  - `P`: A ping-pong matrix that denotes where the latest results are.

  - `rsd`: The `RadixSortDetails` object used to configure the algorithm.
"""
radixsort_par_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, P::AbstractVector{Bool}, rsd::RadixSortDetails) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}} = @inbounds @views begin
  lo = low(rsd)
  hi = high(rsd)

  C = alloc!(rsd, Int(bitmask(rsd))+1)
  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, CountSortDetails(rsd))
  P[lo:hi] .= !P[lo]

  nC = length(C)
  @sync for j in 1:nC
    nlo = C[j]+1
    nhi = j != nC ? C[j+1] : hi

    nrsd = next(rsd, nlo, nhi)

    (is_small(nrsd) || is_deep(nrsd)) && continue

    if is_big(nrsd)
      Threads.@spawn radixsort_par_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nrsd)
    else
      Threads.@spawn radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nrsd)
    end
  end

  free!(rsd, C)
end


"""
$(SIGNATURES)

Parallel-to-Parallel radixsort.

# Arguments
  - `Va`: The cloud of points.
  - `Ra`: The morton transformation for each point of `Va`
  - `Ia`: The current permutation of `Va`.

  - `Vb`: An auxiliary matrix to store the partially sorted version of `Va`.
  - `Rb`: An auxiliary vector to store the partially sorted version of `Ra`.
  - `Ib`: An auxiliary vector to store the resulting permutation.

  - `P`: A ping-pong matrix that denotes where the latest results are.

  - `rsd`: The `RadixSortDetails` object used to configure the algorithm.
"""
radixsort_par_par_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, P::AbstractVector{Bool}, rsd::RadixSortDetails) where {TV<:AbstractMatrix, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}} = @inbounds @views begin
  lo = low(rsd)
  hi = high(rsd)

  Cm = alloc!(rsd, Int(bitmask(rsd))+1, Threads.nthreads())
  C  = alloc!(rsd, Int(bitmask(rsd))+1)
  countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, Cm, CountSortDetails(rsd))
  minimum!(reshape(C, :, 1), Cm)
  free!(rsd, Cm)

  P[lo:hi] .= !P[lo]

  nC = length(C)
  @sync for j in 1:nC
    nlo = C[j]+1
    nhi = j != nC ? C[j+1] : hi
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


