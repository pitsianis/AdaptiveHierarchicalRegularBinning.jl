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

RadixSortDetails(len, par_th, seq_th, sml_th, lvl_th) = RadixSortDetails(1, len, 1, par_th, seq_th, sml_th, lvl_th, [Dict() for _ in 1:Threads.nthreads()])
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


alloc!(a::RadixSortDetails, dims...) = @inbounds begin
  pool = a.pools[Threads.threadid()]
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

free!(a::RadixSortDetails, X::Array{UInt}) = @inbounds begin
  pool = a.pools[Threads.threadid()]
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

radixsort_seq_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd::RadixSortDetails) = @inbounds @views begin

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


radixsort_par_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd::RadixSortDetails) = @inbounds @views begin

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

radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd::RadixSortDetails) = @inbounds @views begin

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


