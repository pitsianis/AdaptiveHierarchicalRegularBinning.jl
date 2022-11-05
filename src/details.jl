const DEFAULT_THRESHOLDS = (PAR = 10_000, SEQ = 1_000, SML=100, DPT=typemax(UInt))

struct CountSortDetails{B, D}
  lo::UInt
  hi::UInt

  l::UInt
end

struct RadixSortDetails{B, D}
  csd::CountSortDetails{B, D}

  par_th::UInt
  seq_th::UInt
  sml_th::UInt
  dpt_th::UInt

  pools::Vector{Dict{Integer, Vector{Vector{UInt}}}}
end


CountSortDetails(bitlen, lo, hi; dims=1) = CountSortDetails{bitlen, dims}(lo, hi, 1)
CountSortDetails(rsd::RadixSortDetails)  = rsd.csd

next(csd::CountSortDetails, lo, hi) = typeof(csd)(lo, hi, depth(csd)+1)


low(csd::CountSortDetails)    = csd.lo
high(csd::CountSortDetails)   = csd.hi
length(csd::CountSortDetails) = high(csd)-low(csd)+1
depth(csd::CountSortDetails)  = csd.l

bitlen(::CountSortDetails{B})         where {B} = B
leaddim(::CountSortDetails{<:Any,D})  where {D} = D

bitmask(csd::CountSortDetails) = (one(UInt) << bitlen(csd)) - 1
radixsel(csd::CountSortDetails, x) = (x >> (8*sizeof(x) - depth(csd)*bitlen(csd))) & bitmask(csd)

@inline function Base.selectdim(A, csd::CountSortDetails, i)
  I = ntuple(k->k==leaddim(csd) ? i : (:), Val(max(ndims(A), leaddim(csd))))
  @boundscheck checkbounds(A, I...)
  return @inbounds @view A[I...]
end


RadixSortDetails(bitlen, lo, hi; dims=1, par_th=DEFAULT_THRESHOLDS.PAR, seq_th=DEFAULT_THRESHOLDS.SEQ, sml_th=DEFAULT_THRESHOLDS.SML, dpt_th=DEFAULT_THRESHOLDS.DPT) = RadixSortDetails{bitlen, dims}(CountSortDetails(bitlen, lo, hi; dims=dims), par_th, seq_th, sml_th, dpt_th, [Dict() for _ in 1:Threads.nthreads()])

next(rsd::RadixSortDetails, lo, hi) = typeof(rsd)(next(rsd.csd, lo, hi), rsd.par_th, rsd.seq_th, rsd.sml_th, rsd.dpt_th, rsd.pools)

for fn in (:low, :high, :length, :depth, :bitlen, :leaddim, :bitmask)
  @eval $fn(rsd::RadixSortDetails) = $fn(CountSortDetails(rsd))
end


is_huge(rsd::RadixSortDetails)  = length(rsd) >= rsd.par_th
is_big(rsd::RadixSortDetails)   = length(rsd) >= rsd.seq_th
is_small(rsd::RadixSortDetails) = length(rsd) <= rsd.sml_th
is_deep(rsd::RadixSortDetails)  = depth(rsd)  >= rsd.dpt_th


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
