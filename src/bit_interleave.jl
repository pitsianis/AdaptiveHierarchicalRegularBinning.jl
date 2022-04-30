using BitOperations

function has_pdep()
  # NOTE: https://discourse.julialang.org/t/bit-manipulation-instruction-set/5368/9
  CPUInfo = zeros(Int32, 4)
  ccall(:jl_cpuidex, Cvoid, (Ptr{Cint}, Cint, Cint), CPUInfo, 7, 0)
  return CPUInfo[2] & 0x100 != 0
end

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
bit_space(::Type{Pdep}, w::UInt8, n::Integer) = convert(UInt8, bit_space(Pdep, convert(UInt32, w), n) & 0xFF)
bit_space(::Type{Pdep}, w::UInt16, n::Integer) = convert(UInt16, bit_space(Pdep, convert(UInt32, w), n) & 0xFFFF)

# Imlementations for Pdep
bit_space(::Type{Pdep}, w::UInt32, n::Integer) = pdep(w, bit_space_mask(typeof(w), n))
bit_space(::Type{Pdep}, w::UInt64, n::Integer) = pdep(w, bit_space_mask(typeof(w), n))

# Generic fallback (mainly for UInt128)
bit_space(::Type{Pdep}, w::Unsigned, n::Integer) = bit_space(Brute, w, n)
# !SECTION
# !SECTION


# SECTION: bit_interleave
# TODO: Document
# TODO: Handle iterators
# TODO: Memoization
function bit_interleave(M::Type{<:InterleaveMethod}, W::AbstractVector{<:Unsigned})
  L = length(W)

  mapper((i, w)) = bit_space(M, w, L - 1) << (i - 1)
  reducer = |

  return mapreduce(mapper, reducer, enumerate(W); init=zero(eltype(W)))
end


bit_interleave(W::AbstractVector{<:Unsigned}) = bit_interleave(@static(has_pdep() ? Pdep : Brute), W)
bit_interleave(Ws::T...) where {T<:Unsigned} = bit_interleave(collect(Ws))
bit_interleave(w::Unsigned) = w

bit_interleave(W::AbstractMatrix{<:Unsigned}; dims::Integer=2) = map(bit_interleave, eachslice(W; dims=dims))
# !SECTION

# SECTION: Helpers
# NOTE: https://lemire.me/blog/2018/01/08/how-fast-can-you-bit-interleave-32-bit-integers/
pdep(x::UInt32, y::UInt32)::UInt32 = ccall("llvm.x86.bmi.pdep.32", llvmcall, UInt32, (UInt32, UInt32), x, y)
pdep(x::UInt64, y::UInt64)::UInt64 = ccall("llvm.x86.bmi.pdep.64", llvmcall, UInt64, (UInt64, UInt64), x, y)

pext(x::UInt32, y::UInt32)::UInt32 = ccall("llvm.x86.bmi.pext.32", llvmcall, UInt32, (UInt32, UInt32), x, y)
pext(x::UInt64, y::UInt64)::UInt64 = ccall("llvm.x86.bmi.pext.64", llvmcall, UInt64, (UInt64, UInt64), x, y)

"""
$(TYPEDSIGNATURES)

Generates a spaced mask.

# Arguments
- `T`: The resulting type
- `n`: The amount of padding

# Examples
```julia-repl
julia> brmask(UInt8, 0) |> bitstring
"11111111"

julia> brmask(UInt8, 1) |> bitstring
"01010101"

julia> brmask(UInt8, 2) |> bitstring
"01001001"

julia> brmask(UInt8, 3) |> bitstring
"00010001"
```
"""
bit_space_mask(T::Type{<:Unsigned}, n::Integer) = bit_space_mask_memoized(T, Val(n))

function bit_space_mask_impl(T::Type{<:Unsigned}, n::Integer)
  # TODO: Pure?
  # TODO: Memoization? https://github.com/JuliaCollections/Memoize.jl
  mask = bmask(T, n)

  # Remaining number of segments
  k = bsizeof(mask) ÷ (n + 1)
  # Shift amount
  s = n + 1

  while k > 0
    mask |= mask << s

    # TODO: Bitshifts?
    s *= 2
    k ÷= 2
  end

  return (mask << 1) | 0x01
end

# FIXME: Memoization hurts performance
@generated bit_space_mask_memoized(::Type{T}, ::Val{n}) where {T<:Unsigned, n} = :($(bit_space_mask_impl(T, n)))
# !SECTION
