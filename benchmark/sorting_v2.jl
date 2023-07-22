using AdaptiveHierarchicalRegularBinning: radixsort_seq_seq_impl!, radixsort_par_seq_impl!, radixsort_par_par_impl!, RadixSortDetails, Allocator, spatial_encode!
using BenchmarkTools


n = 2_000_000
d = 8
l = 8

function sortperm_int_range!(P::Vector{<:Integer}, x::Vector{<:Integer}, rangelen, minval)
  offs = 1 - minval
  n = length(x)

  counts = fill(0, rangelen+1)
  counts[1] = 1
  @inbounds for i = 1:n
      counts[x[i] + offs + 1] += 1
  end

  #cumsum!(counts, counts)
  @inbounds for i = 2:length(counts)
      counts[i] += counts[i-1]
  end

  @inbounds for i = 1:n
      label = x[i] + offs
      P[counts[label]] = i
      counts[label] += 1
  end

  return P
end

function fastpermute!(O::Matrix, I::Matrix, p::Vector)

  @inbounds for j in eachindex(p)
    l = p[j]
    for i in axes(O,1)
      O[i, j] = I[i, l]
    end
  end

end 

function fastpermute!(O::Vector, I::Vector, p::Vector)

  @inbounds for j in eachindex(p)
    l = p[j]
    O[j] = I[l]
  end

end 

function do_benchmark( n, d, l, p )
  V = rand(d, n)
  R = zeros(UInt, n)
  spatial_encode!(R,V,l; dims = Val(2), center = false)

  I = collect(UInt, 1:n)

  Va = copy(V)
  Ra = copy(R)
  Ia = copy(I)

  Vb = similar(V)
  Rb = similar(R)
  Ib = similar(I)

  P = zeros(Bool, n)

  rsd = RadixSortDetails(d, 1, n; dims=2, sml_th=p, dpt_th=l+1)

  allocator = Allocator(eltype(I))

  println("AHRB")

  @btime begin 
    radixsort_seq_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, allocator) 
    
    $Va[:, $P] .= $Vb[:, $P]
    $Ra[$P]    .= $Rb[$P]
    $Ia[$P]    .= $Ib[$P]

  end setup=( Va=deepcopy($Va); Ra=deepcopy($Ra); Ia=deepcopy($Ia); Vb=deepcopy($Vb); Rb=deepcopy($Rb); Ib=deepcopy($Ib); P=deepcopy($P); rsd=deepcopy($rsd); allocator=deepcopy($allocator) )

  println("sortperm")
  # foo = @btime sort($R)

  ix = [LinearIndices(R);];
  rx = copy(R)
  vx = copy(Va)
  
  @btime begin
    # sortperm_int_range!($ix, $Ra, 1024, 0);
    sort!($ix; order = $(Base.Sort.Perm(Base.Order.Forward, vec(R))))
    fastpermute!($rx, $R, $ix)
    fastpermute!($vx, $Va, $ix)
  end

  return nothing

end
  
