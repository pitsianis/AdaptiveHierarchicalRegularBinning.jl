"""
Contains needed information for the radix sort algorithm.

  - `lo`: The lower index.
  - `hi`: The highest index.
  - `lvl`: The level of recursion.

  - `par_th`: Minimum number of elements to use `radixsort_par_par_impl!`
  - `seq_th`: Minimum number of elements to use `radixsort_par_seq_impl!`
  - `sml_th`: Minimum number of elements to sort.
  - `lvl_th`: Maximum level of recursion.

  - `pools`: A memory pool for each thread.
"""
struct RadixSortDetails
  lo::UInt
  hi::UInt
  lvl::UInt

  par_th::UInt
  seq_th::UInt
  sml_th::UInt
  lvl_th::UInt

  pools::Vector{Dict{Integer, Vector{Vector{UInt}}}}
end

"""
$(TYPEDSIGNATURES)

Default initializes a `RadixSortDetails` struct.

  - `len`: The length of the array to sort.
  - `par_th`: The parallel to parallel threshold.
  - `seq_th`: The parallel to sequential threshold.
  - `sml_th`: The small threshold.
  - `lvl_th`: The level threshold.
"""
RadixSortDetails(len::Integer, par_th::Integer, seq_th::Integer, sml_th::Integer, lvl_th::Integer) = RadixSortDetails(1, len, 1, par_th, seq_th, sml_th, lvl_th, [Dict() for _ in 1:Threads.nthreads()])

"""
$(TYPEDSIGNATURES)

Computes the next `RadixSortDetails` object.

  - `rsd`: The previous `RadixSortDetails` object.
  - `nlo`: The next lower boundary.
  - `nhi`: The next highest boundary.
"""
next!(rsd::RadixSortDetails, nlo::UInt, nhi::UInt) = RadixSortDetails(nlo, nhi, level(rsd)+1, rsd.par_th, rsd.seq_th, rsd.sml_th, rsd.lvl_th, rsd.pools)

Base.first(rsd::RadixSortDetails) = rsd.lo
Base.last(rsd::RadixSortDetails)  = rsd.hi
Base.length(rsd::RadixSortDetails) = last(rsd)-first(rsd)+1

level(rsd::RadixSortDetails) = rsd.lvl

is_huge(rsd::RadixSortDetails)  = length(rsd) >= rsd.par_th
is_big(rsd::RadixSortDetails)   = length(rsd) >= rsd.seq_th
is_small(rsd::RadixSortDetails) = length(rsd) <= rsd.sml_th
is_deep(rsd::RadixSortDetails)  = level(rsd)  >= rsd.lvl_th

Base.to_index(rsd::RadixSortDetails) = first(rsd):last(rsd)


"""
$(TYPEDSIGNATURES)

Allocates an `Array` of dimensions `dims`.

  - `rsd`: The `RadixSortDetails` to allocate from.
  - `dims`: The size of the array.
"""
alloc!(rsd::RadixSortDetails, dims...) = @inbounds begin
  pool = rsd.pools[Threads.threadid()]
  len = prod(dims)
  if !haskey(pool, len)
    pool[len] = []
  end

  buffers = pool[len]
  if isempty(buffers)
    push!(buffers, Vector{UInt}(undef, len))
  end

  return reshape(pop!(buffers), dims...)
end

"""
$(TYPEDSIGNATURES)

Frees an `Array` to be reallocated.

  - `rsd`: The `RadixSortDetails` to free to.
  - `X`: The array to free.
"""
free!(rsd::RadixSortDetails, X::Array{UInt}) = @inbounds begin
  pool = rsd.pools[Threads.threadid()]
  dims = size(X)
  len = prod(dims)

  if !haskey(pool, len)
    pool[len] = []
  end

  buffers = pool[len]
  push!(buffers, reshape(X, :))

  return
end



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
  C = alloc!(rsd, 256)

  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, first(rsd), last(rsd), level(rsd))
  P[rsd] .= level(rsd)%2

  for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : last(rsd)
    nrsd = next!(rsd, nlo, nhi)

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

  C = alloc!(rsd, 256)
  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, first(rsd), last(rsd), level(rsd))
  P[rsd] .= level(rsd)%2

  @sync for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : last(rsd)

    nrsd = next!(rsd, nlo, nhi)

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

  Cm = alloc!(rsd, 256, Threads.nthreads())
  C  = alloc!(rsd, 256)
  countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, Cm, first(rsd), last(rsd), level(rsd))
  minimum!(reshape(C, :, 1), Cm)
  free!(rsd, Cm)

  @sync for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : last(rsd)
    nrsd = next!(rsd, nlo, nhi)

    length(nrsd) <= 1 && continue

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


