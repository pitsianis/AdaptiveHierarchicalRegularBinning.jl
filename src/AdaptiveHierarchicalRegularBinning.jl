module AdaptiveHierarchicalRegularBinning

import Base: length, eltype
using DocStringExtensions


include("utilities.jl")
include("bit_interleave.jl")
include("spatial_encode.jl")

include("details.jl")
include("countsort.jl")
include("radixsort.jl")
include("tree.jl")



"""
$(SIGNATURES)

Constructs the tree.

# Arguments
  - `V`: The cloud of points. (`d`x`n`)
  - `l`: The maximum tree depth.
"""
function regural_bin(V, l)
  dims=2 #TODO: Make this arg
  R = Vector{UInt}(undef, size(V, dims))
  bitlen = size(V, dims==1 ? 2 : 1)
  #TODO: Pass information to the tree
  δ, σ = spatial_encode!(R, V, l; dims=Val(dims), center=false)
  I = collect(UInt, 1:length(R))

  Va=V
  Ra=R
  Ia=I

  Vb = similar(V)
  Rb = similar(R)
  Ib = similar(I)

  P = zeros(Bool, length(R))

  rsd = RadixSortDetails(bitlen, 1, length(R); dims=dims, dpt_th=l+1)
  alloc = Allocator(UInt)
  radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, alloc)

  Va[:, P] .= Vb[:, P]
  Ra[P] .= Rb[P]
  Ia[P] .= Ib[P]

  #TODO: Spatial encode should take care of this
  R .= R .<< (sizeof(eltype(R))*8 - bitlen*l)
  tree = make_tree(R, l, bitlen)

  return tree
end

end
