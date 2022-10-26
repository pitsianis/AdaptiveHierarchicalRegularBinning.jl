
# using BitOperations, BenchmarkTools

function translate_scale(x; center = false, dims = 2)
  dims = (dims == 2) ? 1 : 2
  p = copy(x)
  p .-= minimum(p, dims)
  scale = maximum(maximum(p, dims))
  p ./= scale + eps(scale)

  if center
    p .+= minimum(1.0 .- eps(1.0) .- p, dims = 2) ./ 2
  end
  return p
end

"""
$(TYPEDSIGNATURES)

Find the displacement vector δ and scale scalar σ parameters that will
transform the cloud points `p` to be within a unit box [0,1),

0 ≤ σ (p - δ) < 1

- `p`: cloud point coordinates

- `center`: flag to center the cloud by dividing the slack evenly, default false
- `dims`  : dimension of `p` expressing single point coordinates,
i.e. if the `i`^th point is `p[i,:]`, then set `dims=2`.

"""
function translate_scale_vals(p; center = false, dims = 2)

  displ = minimum(p, dims=dims)
  scale = maximum(p, dims=dims)
  scale = maximum(scale .- displ)
  scale += eps(scale)
  scale \= 1.0

  if center
    displ /= 2.0
  end
  return (reshape(displ, :),scale)
end

@inline function quantize(p, L)
  floor.(UInt, p .* 2^L)
end

"""
$(TYPEDSIGNATURES)

Spatially encode a single point `0 ≤ p < 1` for L levels
"""
@inline function morton(p, L)
  t = BitArray((p[i] & 2^j)>0 for i=1:length(p), j = 0:L-1)
  return t.chunks[1]
end

"""
$(TYPEDSIGNATURES)

Spatially encode the point cloud coordinates `x` for L levels
  `dims` denotes the dimension that contains the coordinates of a single point.
"""
function spatial_encode(V, L; dims=2)
  d, n = size(V)
  if dims != 2
    d, n = n, d
  end

  δ, σ = translate_scale_vals(V; dims=dims)

  c = Vector{UInt64}(undef, n)
  @inbounds Threads.@threads for i=1:n
    # FIXME: selectdim -> performance impact
    # Either stick to only one nxd or dxn
    # Or create compile-time selectdim
    x = selectdim(V, dims, i)
    x = reshape(x, :)
    y  = σ * (x .- δ)
    q = quantize.(y, L)
    c[i] = bit_interleave(q)
  end

  return (c, δ, σ)
end


#= ----------- UNDER CONSTRUCTION --------------------
function bit_interleave(v,L)
  c = 0
  d = length(v)
  k = 0
  @inbounds for i = 0:L-1
    for b = 1:d
      c = bset(c,k,bget(v[b],i))
      k += 1
    end
  end
  return c
end

function makeinterleaveLUT()
  T = zeros(UInt32,256,256,256)

  v = zeros(UInt32,3)
  @inbounds for k = 0:255
    v[3] = k
    for j = 0:255
      v[2] = j
      for i = 0:255
        v[1] = i
        T[i+1,j+1,k+1] = bit_interleave(v,8)
      end
    end
  end

  return T
end

@inline interleave8(c1,T) = T[c1+1,1,1]
@inline interleave8(c1,c2,T) = T[c1+1,c2+1,1]
@inline interleave8(c1,c2,c3,T) = T[c1+1,c2+1,c3+1]

function trial(n,L)
  T = makeinterleaveLUT()

  x = rand(3,n)

  @time begin
    mc = zeros(UInt64,n)

    c = quantize(translate_scale(x),L)

    for i = 1:n
      mc[i] =  interleave8(c[:,i]...,T)
    end
  end

  @time begin
    mc2 = zeros(UInt64,n)
    c = quantize(translate_scale(x),L)

    C = BitArray(undef,d,2^L)
    for i = 1:n
      Vector{Char}.(bitstring.(c[i]
      mc2[i] =  interleave8(c[:,i]...,T)
    end

  end

  return
end

function spatial_encode(x,L)
  T = makeinterleaveLUT()

  n = size(x,2)

  @time begin
    displacement, scale = translate_scale_vals(x)
    mc = zeros(UInt64,n)

    c = zeros(UInt32,L)

    for i = 1:n
      c = quantize((x[:,i] .- displacement) .* scale,L)
      mc[i] =  interleave8(c[:]...,T)
    end
  end

  return mc
end

function test(n,d)

  x = randn(d,n)

  @time dsp, scl = translate_scale_vals(x)

  @time c1 = translate_scale(x)
  @time c2 = (x .- dsp) .* scl

  sqrt(sum((c1 .- c2).^2))
end

test(5,4)

x = [0.25 0.25 0.75 0.75; 0.25 0.75 0.25 0.75]
spatial_encode(x,1)

=#
