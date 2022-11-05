struct CountSortDetails{B, D}
  lo::Int
  hi::Int

  l::Int
end


CountSortDetails(bitlen, lo, hi; dims=1) = CountSortDetails{bitlen, dims}(lo, hi, 1)
next(csd::CountSortDetails, lo, hi) = typeof(csd)(lo, hi, depth(csd)+1)


low(csd::CountSortDetails)   = csd.lo
high(csd::CountSortDetails)  = csd.hi
depth(csd::CountSortDetails) = csd.l

bitlen(::CountSortDetails{B})         where {B} = B
leaddim(::CountSortDetails{<:Any,D})  where {D} = D

bitmask(csd::CountSortDetails) = (one(UInt) << bitlen(csd)) - 1
radixsel(csd::CountSortDetails, x) = (x >> (8*sizeof(x) - depth(csd)*bitlen(csd))) & bitmask(csd)

@inline function Base.selectdim(A, csd::CountSortDetails, i)
  I = ntuple(k->k==leaddim(csd) ? i : (:), Val(max(ndims(A), leaddim(csd))))
  @boundscheck checkbounds(A, I...)
  return @inbounds @view A[I...]
end
