using AbstractTrees

#SECTION: Structs

"""
$(TYPEDEF)

Represents a node of the tree.

# Fields
  - `lo`: The start index of the represented group.
  - `hi`: The stop index of the represented group.
  - `dpt`: The depth of the node
  - `nidx`: The index of the node.
  - `pidx`: The index of the parent node.
"""
struct NodeInfo
  lo::Int
  hi::Int
  dpt::Int
  nidx::Int
  pidx::Int
end


"""
$(TYPEDEF)

Represents the tree information.

# Fields
  - `T`: The element type of the points.
  - `E`: The element type of the encoded points.
  - `B`: The bit length of each bit group in the encoded space.
  - `D`: The leading dimension.

  - `points`: The cloud of points that the tree spans.
  - `encoded`: The encoded cloud of `points`.
  - `perm`: The permutation of the `points`.
  - `scale`: The original scale of the dimensions of the `points`.
  - `offset`: The per dimension offset of the `points`.
  - `nodes`: Per node information.
  - `children`: Per node children.
  - `maxdepth`: The maximum depth of the tree.
  - `maxpoints`: Small threshold.
"""
struct TreeInfo{T, E, C, GTC}
  points::Matrix{T}
  encoded::Vector{E}
  perm::Vector{Int}

  scale::T
  offset::Vector{T}

  nodes::Vector{NodeInfo}
  children::Vector{Vector{Int64}}
  context::Vector{C}

  qcenters::Vector{Vector{T}}
  qsidelengths::Vector{T}

  maxdepth::Int
  maxpoints::Int

  globalcontext::GTC
end


"""
$(TYPEDEF)

Represents the tree.

# Fields
  - `T`: The element type of the points.
  - `E`: The element type of the encoded points.
  - `B`: The bit length of each bit group in the encoded space.
  - `D`: The leading dimension.

  - `info`: The tree information.
  - `nidx`: Index to the current root of the tree.
"""
struct SpatialTree{T, E, C, GTC}
  info::TreeInfo{T, E, C, GTC}
  nidx::Int
end
#!SECTION

#SECTION: NodeInfo
NodeInfo(t::SpatialTree) = @inbounds TreeInfo(t).nodes[nindex(t)]

low(n::NodeInfo) = n.lo
high(n::NodeInfo) = n.hi
range(n::NodeInfo) = low(n):high(n)
length(n::NodeInfo) = length(range(n))
depth(n::NodeInfo) = n.dpt
nindex(n::NodeInfo) = n.nidx
pindex(n::NodeInfo) = n.pidx

function Base.show(io::IO, node::NodeInfo)
  print(io, 
  """
  NodeInfo:
  node $(Int(node.nidx)), parent $(Int(node.pidx)), depth $(node.dpt) 
  representing $(Int(node.lo)):$(Int(node.hi)) points""")
end
#!SECTION

#SECTION: TreeInfo
TreeInfo(V, R, I, N, C, depth, maxpoints, scale, offset, centers, sidelengths; ctxtype = Bool, gtctype = Nothing) = TreeInfo{eltype(V), eltype(R), ctxtype, gtctype}(V, R, I, scale, offset, N, C, ctxtype[], centers, sidelengths, depth, maxpoints, gtctype())
TreeInfo(t::SpatialTree) = t.info

eltype(  ::TreeInfo{T}                          ) where T = T
enctype( ::TreeInfo{T, E}       where {T}       ) where E = E
ctxtype( ::TreeInfo{T, E, C}    where {T, E}    ) where C = C
gtctype( ::TreeInfo{T, E, C, B} where {T, E, C} ) where B = B
#!SECTION

#SECTION: SpatialTree
SpatialTree(info::TreeInfo, nidx) = SpatialTree{eltype(info), enctype(info), ctxtype(info), gtctype(info)}(info, nidx)

nindex(t::SpatialTree) = t.nidx
cindices(t::SpatialTree) = @inbounds TreeInfo(t).children[nindex(t)]

for fn in (:range, :low, :high, :length, :depth, :pindex)
  @eval $fn(t::SpatialTree) = $fn(NodeInfo(t))
end

for fn in (:enctype, :eltype, :ctxtype, :gtctype)
  @eval $fn(t::SpatialTree) = $fn(TreeInfo(t))
end


points(t::SpatialTree) = @inbounds view(TreeInfo(t).points, :, range(t))
encpoints(t::SpatialTree) = @inbounds @view TreeInfo(t).encoded[range(t)]
isdeep(t::SpatialTree) = depth(t) >= TreeInfo(t).maxdepth
issmall(t::SpatialTree) = length(t) <= TreeInfo(t).maxpoints
qcenter(t::SpatialTree) = @inbounds TreeInfo(t).qcenters[nindex(t)]
qsidelength(t::SpatialTree) = @inbounds TreeInfo(t).qsidelengths[nindex(t)]
const isleaf = isempty ∘ cindices

"""
    leafcount(node)

Get the number of leaves of the tree rooted at `node`.

"""
# this recurses through all nodes in the tree and so may be slow
leafcount(node) = isleaf(node) ? 1 : mapreduce(leafcount, +, children(node))

# user-provided context per node
"""
$(TYPEDEF)

Set context of a node
"""
function setcontext!(t::SpatialTree{T,E,C,GTC}, v::C) where {T,E,C,GTC}
  TreeInfo(t).context[nindex(NodeInfo(t))] = v
end

"""
$(TYPEDEF)

Read context of a node
"""
getcontext(t::SpatialTree) = TreeInfo(t).context[nindex(NodeInfo(t))]

getglobalcontext(t::SpatialTree) = TreeInfo(t).globalcontext

function Base.show(io::IO, tree::SpatialTree)
  print(io, 
  """
  SpatialTree: 
  $(typeof(tree.info.points))($(size(points(tree))[1]),$(Int(size(points(tree))[2]))) points
  $(treesize(tree)) nodes, $(leafcount(tree)) leaves and max depth $(treeheight(tree))""")
end
#!SECTION

#SECTION: AbstractTrees
AbstractTrees.getroot(t::SpatialTree)   = SpatialTree(TreeInfo(t), 1)
AbstractTrees.parent(t::SpatialTree)    = SpatialTree(TreeInfo(t), pindex(t))
function AbstractTrees.children(t::SpatialTree; by=nothing)
  children = map((i)->SpatialTree(TreeInfo(t), i), cindices(t))
  by === nothing && return children
  sort!(children; by=by)
end
AbstractTrees.nodevalue(t::SpatialTree) = first(encpoints(t))
AbstractTrees.isroot(t::SpatialTree)    = nindex(t) == 1
#!SECTION


scalar(tree::SpatialTree) = tree.info.scale
translation(tree::SpatialTree) = tree.info.offset

_qcenter!(center, node::SpatialTree) = _qcenter!(center, nodevalue(node), size(points(node),1), depth(node), node.info.maxdepth)
_qcenter(T, node::SpatialTree) = _qcenter!(Vector{T}( undef, size(points(node),1) ), node)
_qcenter(node::SpatialTree) = _qcenter(eltype(node), node)


function _qcenter!(center, tag, bitlen, depth, maxdepth)
  fill!(center, 0.5)
  @inbounds for k = 1:maxdepth
    l = maxdepth-k+1
    o = 2.0^(-l-1)
    for bit in 1:bitlen
      if l <= depth
        center[bit] += ifelse(tag & 1 != 0, o, -o)
      end
      tag >>= 1
    end 
  end

  return center
end

function center!(center, node::SpatialTree)
  center .= qcenter(node) .* scalar(node) .+ translation(node)
  return center
end
center(T, node::SpatialTree) = center!(Vector{T}(undef, size(points(node),1)), node)
center(node::SpatialTree) = center(eltype(node), node)


_qsidelength(T, node::SpatialTree)  = _qsidelength(T, depth(node))
_qsidelength(node::SpatialTree) = _qsidelength(eltype(node), depth(node))

_qsidelength(T, depth) = one(T)/(1<<(depth))
_qsidelength(depth) = _qsidelength(typeof(eps()), depth)


sidelength(T, node::SpatialTree) = T(qsidelength(node) * scalar(node))
sidelength(node::SpatialTree) = sidelength(eltype(node), node)

coincidence(node::SpatialTree, other::SpatialTree) = depth(node) == depth(other) && nodevalue(node) == nodevalue(other)
