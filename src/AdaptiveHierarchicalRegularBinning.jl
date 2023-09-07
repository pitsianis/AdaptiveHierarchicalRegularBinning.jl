module AdaptiveHierarchicalRegularBinning

import Base: length, eltype, range
using DocStringExtensions
using Transducers, LoopVectorization, ThreadsX
using DataStructures
using SparseArrays

# TODO: Export less things
export SpatialTree, TreeInfo, NodeInfo, ahrb, nindex, cindices
export range, low, high, depth, pindex, bitlen, enctype, leaddim, eltype, isleaf, leafcount
export points, encpoints, isdeep, qcenter, center, qsidelength, sidelength, staticselectdim
export original_perm, original_perm!
export setcontext!, getcontext, getglobalcontext
export qbox2boxdistInf, qbox2boxdist, box2boxdist, point2boxdist
export adjacentLeaf2LeafAnyLevel, adjacentSameLevel, neighborhood
export dualtreetraversal, multilevelinteractions, prioritymultilevelinteractions

include("tree_structures.jl")

include("utilities.jl")

include("spatial_encode.jl")
include("count_and_permutation.jl")
include("build_tree.jl")

include("boxdistance.jl")
include("dualtreetraversal.jl")

include("blockecp.jl")
include("fixedlengthencoding.jl")

include("knn.jl")

function __init__()

  # Decide whether to enable multithreading or not
  if Threads.nthreads() > 1
    enable_multithreading()
  else
    disable_multithreading()
  end

end


"""
    $(TYPEDSIGNATURES)

Generate an adaptive hierarchical regular binning tree structure.

# Arguments
- `V`: A matrix of data points, where each column represents a data point and each row represents a
  feature.
- `maxdepth`: An integer that specifies the maximum depth of the tree. Default is
  `ceil(Int, log2(size(V,2)))`.
- `maxpoints`: An integer that specifies the maximum number of points in a leaf node. Default is
  `ceil(Int, (size(V,2))^(1/4))`.
- `QT`: (Optional) The type of the vector `R` used to store the spatial encoding of the data points.
  Defaults to `UInt`.
- `ctxtype`: (Optional) The type of the context information stored at each node of the tree.
  Defaults to `Nothing`.
- `gtctype`: (Optional) The type of the global context information stored alongside the tree.
  Defaults to `Nothing`.
- `lstep`: (Optional) The number of levels to process in each block iteration. Defaults to `2`.
- `method`: (Optional) The method to use for the spatial encoding. Defaults to `"block-ecp"`.
  Possible values are:
  - `"block-ecp"`: Use the block ECP version.
  - `"fixed-length"`: Use the fixed length encoding.

# Returns
The generated tree structure.

# Example
```jldoctest; filter = [r".*nodes,.*leaves and max depth.*"]
julia> V = rand(6, 10_000);

julia> tree = ahrb(V, 5, 10)
SpatialTree: 
Matrix{Float64}(6,10000) points
3805 nodes, 3740 leaves and max depth 2
```
"""
function ahrb(V, maxdepth = ceil(Int, log2(size(V,2))), maxpoints = ceil(Int, (size(V,2))^(1/4)); QT = UInt, ctxtype = Nothing, gtctype = Nothing, method = "block-ecp", lstep::Int = 2)::SpatialTree
  if method == "block-ecp"
    return ahrb_block(V, maxdepth, maxpoints; QT, ctxtype, gtctype, lstep)
  elseif method == "fixed-length"
    return ahrb_fixed_length(V, maxdepth, maxpoints; QT, ctxtype, gtctype)
  else
    error("Unknown method: $method")
  end
end

end
