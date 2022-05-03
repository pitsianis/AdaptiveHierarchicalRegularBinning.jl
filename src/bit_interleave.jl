using BitOperations

# SECTION: Precompilation dependencies
function has_pdep()
  # NOTE: https://discourse.julialang.org/t/bit-manipulation-instruction-set/5368/9
  CPUInfo = zeros(Int32, 4)
  ccall(:jl_cpuidex, Cvoid, (Ptr{Cint}, Cint, Cint), CPUInfo, 7, 0)
  return CPUInfo[2] & 0x100 != 0
end
# !SECTION

# SECTION: Method Selectors
abstract type InterleaveMethod end
abstract type Brute <: InterleaveMethod end
abstract type Pdep <: InterleaveMethod end
# !SECTION

# SECTION: bit_space
"""
$(FUNCTIONNAME)(M, w, n)

Evenly space out bits.

# Arguments
- `M`: The method to space out the bits
- `w`: The word to space out
- `n`: The amount of padding

# Examples
```julia-repl
julia> bit_space(Brute, 0x00ABCDEF, 0) |> bitstring
"00000000101010111100110111101111"

julia> bit_space(Brute, 0x00ABCDEF, 1) |> bitstring
"10100000101000101010100010101010"

julia> bit_space(Brute, 0x00ABCDEF, 2) |> bitstring
"00000100100100100000100100100100"

julia> bit_space(Brute, 0x00ABCDEF, 3) |> bitstring
"10001000100000001000100010001000"

julia> bit_space(Pdep, 0x00ABCDEF, 2) == bit_space(Brute, 0x00ABCDEF, 2)
true
```
"""
bit_space(M::Type{<:InterleaveMethod}, w::Unsigned, n::Integer) = throw("$bit_space is not implemented for method $M")

# SECTION: Brute
function bit_space(::Type{Brute}, w::Unsigned, n::Integer)
  # TODO: Dynamic type expansion?
  #       bsizeof(ret) >= L*bsizeof(w)
  mapper(i) = ((w >> i) & 0x1) << (i * (n + 1))

  return mapreduce(mapper, |, 0:bsizeof(w)-1; init=zero(eltype(w)))
end
# !SECTION

# SECTION: pdep
# Conversions for small UInts
bit_space(::Type{Pdep}, w::UInt8,  n::Integer) = convert(UInt8,  bit_space(Pdep, convert(UInt32, w), n) & 0xFF)
bit_space(::Type{Pdep}, w::UInt16, n::Integer) = convert(UInt16, bit_space(Pdep, convert(UInt32, w), n) & 0xFFFF)

# Imlementations for Pdep
bit_space(::Type{Pdep}, w::UInt32, n::Integer) = pdep(w, bit_space_mask(UInt32, n))
bit_space(::Type{Pdep}, w::UInt64, n::Integer) = pdep(w, bit_space_mask(UInt64, n))

# Generic fallback (mainly for UInt128)
bit_space(::Type{Pdep}, w::Unsigned, n::Integer) = bit_space(Brute, w, n)
# !SECTION
# !SECTION


# SECTION: bit_interleave
# TODO: Handle known-sized iterators
# TODO: Memoization
"""
    $(FUNCTIONNAME)([M, ]W)

Interleaves an array of unsigned integers.

# Arguments
- `M`: The method with which to interleave
- `W`: The array

# Examples
```julia-repl
julia> $(FUNCTIONNAME)([0x00FF, 0x000F]) |> bitstring
"0101010111111111"

julia> $(FUNCTIONNAME)([0x000F, 0x00FF]) |> bitstring
"1010101011111111"

julia> $(FUNCTIONNAME)([0x0080, 0x0001]) |> bitstring
"0100000000000010"

julia> $(FUNCTIONNAME)([0x0001, 0x0080]) |> bitstring
"1000000000000001"
```
"""
function bit_interleave(M::Type{<:InterleaveMethod}, W::AbstractVector{<:Unsigned})
  # TODO: W does not have to be AbstractVector (Any iterator with constant eltype and known size is ok)
  L = length(W)

  L == 0 && return zero(eltype(W))
  L == 1 && return first(W)

  # TODO: mapreduce is slower than this implementation
  ret = zero(eltype(W))
  for w in Iterators.reverse(W)
    ret = (ret << 1) | bit_space(M, w, L-1)
  end

  return ret
end

bit_interleave(W::AbstractVector{<:Unsigned}) = bit_interleave(@static(has_pdep() ? Pdep : Brute), W)

"""
$(SIGNATURES)

Interleaves all unsigned integers.
"""
bit_interleave(Ws::T...) where {T<:Unsigned} = bit_interleave(collect(Ws))

"""
    $(FUNCTIONNAME)(W[; dims])

Interleaves all unsigned integers along a dimention.

# Arguments
- `W`: The matrix
- `dims`: The dimention. Defaults to 1
"""
bit_interleave(W::AbstractMatrix{<:Unsigned}; dims::Integer=1) = map(bit_interleave, eachslice(W; dims=dims))
# !SECTION

# SECTION: Helpers
# NOTE: https://lemire.me/blog/2018/01/08/how-fast-can-you-bit-interleave-32-bit-integers/
# SECTION: pdep
pdep(x::UInt32, y::UInt32)::UInt32 = ccall("llvm.x86.bmi.pdep.32", llvmcall, UInt32, (UInt32, UInt32), x, y)
pdep(x::UInt64, y::UInt64)::UInt64 = ccall("llvm.x86.bmi.pdep.64", llvmcall, UInt64, (UInt64, UInt64), x, y)

pext(x::UInt32, y::UInt32)::UInt32 = ccall("llvm.x86.bmi.pext.32", llvmcall, UInt32, (UInt32, UInt32), x, y)
pext(x::UInt64, y::UInt64)::UInt64 = ccall("llvm.x86.bmi.pext.64", llvmcall, UInt64, (UInt64, UInt64), x, y)
# !SECTION


# SECTION: bit_space_mask
function bit_space_mask_impl(T::Type{<:Unsigned}, n::Integer)
  mask = bmask(T, n)
  # Shift amount
  s = n+1

  # Last segment mask
  lsmask = bmask(T, bsizeof(T)-s:bsizeof(T))
  lsmask == 0 && return one(T)

  while mask & lsmask == 0
    mask |= mask << s
    s <<= 1
  end

  return (mask << 1) | 0x01
end

"""
$(SIGNATURES)

Generates a spaced mask.

# Arguments
- `T`: The resulting type
- `n`: The amount of padding

# Examples
```julia-repl
julia> $(FUNCTIONNAME)(UInt8, 0) |> bitstring
"11111111"

julia> $(FUNCTIONNAME)(UInt8, 1) |> bitstring
"01010101"

julia> $(FUNCTIONNAME)(UInt8, 2) |> bitstring
"01001001"

julia> $(FUNCTIONNAME)(UInt8, 3) |> bitstring
"00010001"
```
"""
bit_space_mask(::Type{<:Unsigned}, ::Integer) = throw("Not Implemented")

for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
  let
    # Memoization
    cache = ntuple((i)->bit_space_mask_impl(T, i-1), bsizeof(T))
    eval(:(bit_space_mask(::Type{$T}, n::Integer) = n >= $(length(cache)) ? $(last(cache)) : @inbounds $(cache)[n+1]))
  end
end
# !SECTION
# !SECTION
