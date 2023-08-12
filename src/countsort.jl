"""
$(SIGNATURES)

Sequential countsort.

# Remarks
Inputs `Va`, `Ra`, `Ia` are not altered.

# Arguments
  - `Va`: The cloud of points.
  - `Ra`: The morton transformation for each point of `Va`
  - `Ia`: The current permutation of `Va`.

  - `Vb`: A matrix to store the sorted version of `Va`.
  - `Rb`: A vector to store the sorted version of `Ra`.
  - `Ib`: A vector to store the resulting permutation.

  - `C`: A vector acting as a working space for the countsort algorithm.

  - `csd`: The `CountSortDetails` object used to configure the algorithm.
"""
countsort_seq_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC, csd::CountSortDetails) where {TV<:AbstractArray, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractVector{<:Unsigned}} = @inbounds @views begin
  lo = low(csd)
  hi = high(csd)

  to_index(x) = radixsel(csd, x) + 1

  fill!(C, 0)
  for i in lo:hi
    C[to_index(Ra[i])] += 1
  end

  C[1] += lo - 1
  cumsum!(C, C)

  for i in lo:hi
    ci = to_index(Ra[i])

    j = C[ci]
    C[ci] -= 1

    selectdim(Vb, csd, j) .= selectdim(Va, csd, i)
    Rb[j] = Ra[i]
    Ib[j] = Ia[i]
  end
end


"""
$(SIGNATURES)

Parallel countsort.

# Remarks
Inputs `Va`, `Ra`, `Ia` are not altered.

# Arguments
  - `Va`: The cloud of points.
  - `Ra`: The morton transformation for each point of `Va`
  - `Ia`: The current permutation of `Va`.

  - `Vb`: A matrix to store the sorted version of `Va`.
  - `Rb`: A vector to store the sorted version of `Ra`.
  - `Ib`: A vector to store the resulting permutation.

  - `C`: A matrix acting as a working space for the countsort algorithm.

  - `csd`: The `CountSortDetails` object used to configure the algorithm.
"""
countsort_par_impl!(Va::TV, Ra::TR, Ia::TI, Vb::TV, Rb::TR, Ib::TI, C::TC, csd::CountSortDetails) where {TV<:AbstractArray, TR<:AbstractVector{<:Unsigned}, TI<:AbstractVector{<:Unsigned}, TC<:AbstractMatrix{<:Unsigned}} = @inbounds @views begin
  lo = low(csd)
  hi = high(csd)

  n  = length(csd)
  np = Threads.nthreads()
  b  = cld(n, np)

  to_index(x) = radixsel(csd, x) + 1
  local_range(k) = (k-1)*b + lo : min(k*b + lo - 1, hi)

  Threads.@threads for k=1:np
    fill!(C[:, k], 0)
    for i in local_range(k)
      C[to_index(Ra[i]), k] += 1
    end
  end

  C[1] += lo - 1
  cumsum!(C, C, dims=2)
  cumsum!(C[:, end], C[:, end])
  C[2:end, 1:end-1] .+= C[1:end-1, end]

  Threads.@threads for k=1:np
    for i in local_range(k)
      ci = to_index(Ra[i])

      j = C[ci, k]
      C[ci, k] -= 1

      selectdim(Vb, csd, j) .= selectdim(Va, csd, i)
      Rb[j] = Ra[i]
      Ib[j] = Ia[i]
    end
  end
end
