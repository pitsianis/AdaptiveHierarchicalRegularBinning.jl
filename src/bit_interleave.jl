using BitOperations

# SECTION: Precompilation dependencies
function has_pdep()
  try
    # NOTE: https://discourse.julialang.org/t/bit-manipulation-instruction-set/5368/9
    CPUInfo = zeros(Int32, 4)
    ccall(:jl_cpuidex, Cvoid, (Ptr{Cint}, Cint, Cint), CPUInfo, 7, 0)
    return CPUInfo[2] & 0x100 != 0
  catch
    return false
  end
end
# !SECTION

# SECTION: bit_space
# SECTION: Brute
function bit_space_brute(w::Unsigned, n::Integer)
  # TODO: Dynamic type expansion?
  #       bsizeof(ret) >= L*bsizeof(w)
  mapper(i) = ((w >> i) & 0x1) << (i * (n + 1))

  return mapreduce(mapper, |, 0:bsizeof(w)-1; init=zero(eltype(w)))
end
# !SECTION

# SECTION: pdep
# Conversions for small UInts
bit_space_pdep(w::T, n::Integer) where {T <: Union{UInt8, UInt16}} =
    convert(UInt32, w) |>
    (x) -> bit_space_pdep(x, n) |>
    (x) -> convert(T, x & typemax(T))

# Implementations for Pdep
bit_space_pdep(w::T, n::Integer) where {T <: Union{UInt32, UInt64}} = pdep(w, bit_space_mask(T, n))

# Generic fallback (mainly for UInt128)
bit_space_pdep(w::Unsigned, n::Integer) = bit_space_brute(w, n)
# !SECTION

"""
$(FUNCTIONNAME)(w, n)

Evenly space out bits.

# Arguments
- `w`: The word to space out
- `n`: The amount of padding

# Examples
```julia-repl
julia> $(FUNCTIONNAME)(0x00ABCDEF, 0) |> bitstring
"00000000101010111100110111101111"

julia> $(FUNCTIONNAME)(0x00ABCDEF, 1) |> bitstring
"10100000101000101010100010101010"

julia> $(FUNCTIONNAME)(0x00ABCDEF, 2) |> bitstring
"00000100100100100000100100100100"

julia> $(FUNCTIONNAME)(0x00ABCDEF, 3) |> bitstring
"10001000100000001000100010001000"
```
"""
bit_space(w::Unsigned, n::Integer) = @static(has_pdep() ? bit_space_pdep(w, n) : bit_space_brute(w, n))
# !SECTION


# SECTION: bit_interleave
# TODO: Handle known-sized iterators
# TODO: Memoization
"""
    $(FUNCTIONNAME)(W)

Interleaves an array of unsigned integers.

# Arguments
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
function bit_interleave(W::AbstractVector{<:Unsigned})
  # TODO: W does not have to be AbstractVector (Any iterator with constant eltype and known size is ok)
  L = length(W)

  L == 0 && return zero(eltype(W))
  L == 1 && return first(W)

  # TODO: mapreduce is slower than this implementation
  ret = zero(eltype(W))
  for w in Iterators.reverse(W)
    ret = (ret << 1) | bit_space(w, L-1)
  end

  return ret
end

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
