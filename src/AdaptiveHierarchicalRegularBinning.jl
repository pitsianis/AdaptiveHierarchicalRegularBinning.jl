module AdaptiveHierarchicalRegularBinning

import Base: length, eltype, range
using DocStringExtensions

# TODO: Export less things
export SpatialTree, TreeInfo, NodeInfo, ahrb!, nindex, cindices
export range, low, high, depth, pindex, bitlen, enctype, leaddim, eltype, isleaf
export points, encpoints, isdeep, qcenter, center, qsidelength, sidelength, staticselectdim
export original_perm, original_perm!
export applypostorder!, applypreorder!
export setcontext!, getcontext
export coincidence

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
  - `V`: The cloud of points.
  - `maxL`: The maximum tree depth.
  - `maxP`: Small threshold.

# Keyword Arguments
  - `dims`: dimension than enumerates the points. Defaults to 2.
  - `QT`: datype for the quantized coordinate word vector. Defaults to UInt
"""
function ahrb!(V, maxL, maxP; dims = 2, QT = UInt)
  R = Vector{QT}(undef, size(V, dims))
  bitlen = size(V, dims==1 ? 2 : 1)

  maxL * bitlen > sizeof(QT) * 8 && throw("Not enough bits to represent the requested tree")

  offset, scale = spatial_encode!(R, V, maxL; dims=Val(dims), center=false)
  #TODO: Spatial encode should take care of this
  R .= R .<< (sizeof(eltype(R))*8 - bitlen*maxL)
  I = collect(UInt, 1:length(R))

  Va=V
  Ra=R
  Ia=I

  Vb = similar(V)
  Rb = similar(R)
  Ib = similar(I)

  P = zeros(Bool, length(R))

  # TODO: Have this as an argument
  # Constructs the tree with a 16bit radix to save on memory.
  rbitlen = 16
  rdpt    = cld(bitlen*maxL, rbitlen)

  # TODO: Pass thresholds as parameters
  rsd = RadixSortDetails(rbitlen, 1, length(R); dims=dims, dpt_th=rdpt, sml_th=maxP)
  alloc = Allocator(UInt)
  radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, alloc)

  selectdim(Va, dims, P) .= selectdim(Vb, dims, P)
  Ra[P] .= Rb[P]
  Ia[P] .= Ib[P]

  tree = make_tree(V, R, I, maxL, maxP, bitlen, scale, offset; dims=dims)

  return tree
end

end
