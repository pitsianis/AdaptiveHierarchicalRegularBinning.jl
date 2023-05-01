const DEFAULT_THRESHOLDS = (PAR=100_000, SEQ=10_000, SML=100, DPT=typemax(Int))

"""
$(TYPEDEF)

Contains all of the details needed for the `countsort_xxx_impl!` algorithms.

# Fields
  - `B`: The length of the bit mask to use.
  - `D`: The leading dimension.
  - `lo`: Lower array bound.
  - `hi`: Higher array bound.
  - `l`: The level of recursion.
"""
struct CountSortDetails{B, D}
  lo::Int
  hi::Int

  l::Int
end


"""
$(TYPEDEF)

Contains all of the details needed for the `radixsort_xxx_xxx_impl!` algorithms.

# Fields
  - `B`: The length of the bit mask to use.
  - `D`: The leading dimension.
  - `csd`: The underlying `CountSortDetails` object.

  - `par_th`: Threshold for parallel-to-parallel execution.
  - `seq_th`: Threshold for parallel-to-sequential execution.
  - `sml_th`: Threshold for skiping small arrays.
  - `dpt_th`: recursion threshold.
"""
struct RadixSortDetails{B, D}
  csd::CountSortDetails{B, D}

  par_th::Int
  seq_th::Int
  sml_th::Int
  dpt_th::Int
end


# SECTION: CountSortDetails
"""
$(SIGNATURES)

Constructs a new `CountSortDetails` object.

# Arguments
  - `bitlen`: The bit length of the bit mask used by the algorithm.
  - `lo`: Lower array bound.
  - `hi`: Higher array bound.

# Keyword Arguments
  - `dims`: The leading dimension.
"""
CountSortDetails(bitlen, lo, hi; dims=1) = CountSortDetails{bitlen, dims}(lo, hi, 1)

"""
$(SIGNATURES)

Gets a `CountSortDetails` object from a `RadixSortDetails` object.
"""
CountSortDetails(rsd::RadixSortDetails) = rsd.csd

"""
$(SIGNATURES)

Constructs a new `CountSortDetails` object with a deeper level of recursion.

# Arguments
  - `lo`: New lower array bound.
  - `hi`: New higher array bound.
"""
next(csd::CountSortDetails, lo, hi) = typeof(csd)(lo, hi, depth(csd)+1)


"""
$(SIGNATURES)

The index type of the struct.
"""
eltype(csd::CountSortDetails) = typeof(csd.hi)


"""
$(SIGNATURES)

Gets the lower index of the object.
"""
low(csd::CountSortDetails) = csd.lo

"""
$(SIGNATURES)

Gets the higher index of the object.
"""
high(csd::CountSortDetails) = csd.hi

"""
$(SIGNATURES)

Gets the length of the object.
"""
length(csd::CountSortDetails) = high(csd)-low(csd)+1

"""
$(SIGNATURES)

Gets the level of recursion of the object.
"""
depth(csd::CountSortDetails) = csd.l


"""
$(SIGNATURES)

Gets the bit length the object.
"""
bitlen(::CountSortDetails{B}) where {B} = B

"""
$(SIGNATURES)

Gets the leading dimension the object.
"""
leaddim(::CountSortDetails{<:Any,D}) where {D} = D

"""
$(SIGNATURES)

Gets the bit mask the object.
"""
bitmask(csd::CountSortDetails) = (one(UInt) << bitlen(csd)) - 1

"""
$(SIGNATURES)

Selects the corresponding radix.

# Arguments
  - `csd`: The `CountSortDetails` object to base the selection on.
  - `x`: The value from which the radix is computed.
"""
radixsel(csd::CountSortDetails, x) = radixshft(x, depth(csd), bitlen(csd)) & bitmask(csd)


"""
$(SIGNATURES)

Shifts the corresponding amount of bits.

# Arguments
  - `x`: The value from which the radix is computed.
  - `depth`: The current depth.
  - `bitlen`: The bit length of the bit groups.
"""
radixshft(x, depth, bitlen) = (x >> (8*sizeof(x) - depth*bitlen))

@inline Base.selectdim(A, csd::CountSortDetails, i) = staticselectdim(A, Val(leaddim(csd)), i)
# !SECTION: CountSortDetails

# SECTION: RadixSortDetails
"""
$(SIGNATURES)

Constructs a new `RadixSortDetails` object.

# Arguments
  - `bitlen`: The bit length of the bit mask used by the algorithm.
  - `lo`: Lower array bound.
  - `hi`: Higher array bound.

# Keyword Arguments
  - `par_th`: Threshold for parallel-to-parallel execution.
  - `seq_th`: Threshold for parallel-to-sequential execution.
  - `sml_th`: Threshold for skiping small arrays.
  - `dpt_th`: recursion threshold.
"""
RadixSortDetails(bitlen, lo, hi; dims=1, par_th=DEFAULT_THRESHOLDS.PAR,
                                         seq_th=DEFAULT_THRESHOLDS.SEQ,
                                         sml_th=DEFAULT_THRESHOLDS.SML,
                                         dpt_th=DEFAULT_THRESHOLDS.DPT) =
  RadixSortDetails{bitlen, dims}(CountSortDetails(bitlen, lo, hi; dims=dims), par_th, seq_th, sml_th, dpt_th)

"""
$(SIGNATURES)

Constructs a new `RadixSortDetails` object with a deeper level of recursion.

# Arguments
  - `lo`: New lower array bound.
  - `hi`: New higher array bound.
"""
next(rsd::RadixSortDetails, lo, hi) = typeof(rsd)(next(rsd.csd, lo, hi), rsd.par_th, rsd.seq_th, rsd.sml_th, rsd.dpt_th)

for fn in (:eltype, :low, :high, :length, :depth, :bitlen, :leaddim, :bitmask)
  @eval $fn(rsd::RadixSortDetails) = $fn(CountSortDetails(rsd))
end


is_huge(rsd::RadixSortDetails)  = length(rsd) >= rsd.par_th
is_big(rsd::RadixSortDetails)   = length(rsd) >= rsd.seq_th
is_small(rsd::RadixSortDetails) = length(rsd) <= rsd.sml_th
is_deep(rsd::RadixSortDetails)  = depth(rsd)  >= rsd.dpt_th
# !SECTION: RadixSortDetails
