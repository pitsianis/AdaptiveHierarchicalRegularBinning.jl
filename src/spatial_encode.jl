"""
$(SIGNATURES)

Compute displacement vector and the scale to transform the cloud to a unit hypercube.

# Arguments
  - `V`: The input matrix.

# Keyword Arguments
  - `dims`: The leading dimension.
  - `center`: Whether to divide the slack evenly.
"""
function translate_scale_vals(V::AbstractMatrix; dims, center)
  displ = minimum(V, dims=dims)
  scale = maximum(V, dims=dims)
  scale = maximum(scale .- displ)
  scale += eps(scale)
  scale \= 1.0

  if center
    displ /= 2.0
  end
  return (reshape(displ, :), scale)
end


"""
$(SIGNATURES)

Quantizes a single coordinate.

# Arguments
  - `T`: The resulting type.
  - `x`: The input coordinate.
  - `l`: The levels of the quantizer.
"""
function quantize(T, x, l)
  return floor(T, x * (2^l))
end


"""
$(SIGNATURES)

Encodes a set of points using the morton encoding.

# Arguments
  - `R`: The resulting vector.
  - `V`: The input matrix. Will be scaled and translated.
  - `l`: The levels of the encoding.

# Keyword Arguments
  - `dims`: The leading dimension.
  - `center`: Whether to divide the slack evenly.
"""
function spatial_encode!(R::AbstractVector{<:Integer}, V::AbstractMatrix, l::Integer; dims::Val{D}, center) where {D}
  Di = ndims(V) - (D-1)
  δ, σ = translate_scale_vals(V; dims=D, center=center)

  n = size(V, D)
  nT = Threads.nthreads()
  b = cld(n, nT)

  @inline local_quantizer(x) = quantize(eltype(R), x, l)
  @inline local_range(k) = (k-1)*b + 1 : min(k*b, n)

  Threads.@threads for k in 1:nT
    q = Vector{eltype(R)}(undef, size(V, Di))
    v = Vector{eltype(V)}(undef, size(V, Di))
    @inbounds for i in local_range(k)
      v .= σ .* (staticselectdim(V, dims, i) .- δ)
      q .= local_quantizer.(v)
      R[i] = bit_interleave(q)
    end
  end

  return (δ, σ)
end
