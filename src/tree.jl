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
  lo::UInt
  hi::UInt
  dpt::UInt
  nidx::UInt
  pidx::UInt
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
  - `smlth`: Small threshold.
"""
struct TreeInfo{T, E, B, D}
  points::Matrix{T}
  encoded::Vector{E}
  perm::Vector{UInt}

  scale::T
  offset::Vector{T}

  nodes::Vector{NodeInfo}
  children::Vector{Vector{UInt}}

  maxdepth::UInt
  smlth::UInt
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
struct SpatialTree{T, E, B, D}
  info::TreeInfo{T, E, B, D}
  nidx::UInt
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
#!SECTION

#SECTION: TreeInfo
TreeInfo(V, R, I, N, C, bitlen, depth, smlth, scale, offset; dims) = TreeInfo{eltype(V), eltype(R), bitlen, dims}(V, R, I, scale, offset, N, C, depth, smlth)
TreeInfo(t::SpatialTree) = t.info

eltype(  ::TreeInfo{T}                          ) where T = T
enctype( ::TreeInfo{T, E}       where {T}       ) where E = E
bitlen(  ::TreeInfo{T, E, B}    where {T, E}    ) where B = B
leaddim( ::TreeInfo{T, E, B, D} where {T, E, B} ) where D = D
#!SECTION

#SECTION: SpatialTree
SpatialTree(info::TreeInfo, nidx) = SpatialTree{eltype(info), enctype(info), bitlen(info), leaddim(info)}(info, nidx)

nindex(t::SpatialTree) = t.nidx
cindices(t::SpatialTree) = @inbounds TreeInfo(t).children[nindex(t)]

for fn in (:range, :low, :high, :length, :depth, :pindex)
  @eval $fn(t::SpatialTree) = $fn(NodeInfo(t))
end

for fn in (:bitlen, :enctype, :leaddim, :eltype)
  @eval $fn(t::SpatialTree) = $fn(TreeInfo(t))
end


points(t::SpatialTree) = @inbounds staticselectdim(TreeInfo(t).points, Val(leaddim(t)), range(t))
encpoints(t::SpatialTree) = @inbounds @view TreeInfo(t).encoded[range(t)]
isdeep(t::SpatialTree) = depth(t) >= TreeInfo(t).maxdepth
issmall(t::SpatialTree) = length(t) <= TreeInfo(t).smlth
#!SECTION

#SECTION: AbstractTrees
AbstractTrees.getroot(t::SpatialTree)   = SpatialTree(TreeInfo(t), 1)
AbstractTrees.parent(t::SpatialTree)    = SpatialTree(TreeInfo(t), pindex(t))
AbstractTrees.children(t::SpatialTree)  = Iterators.map((i)->SpatialTree(TreeInfo(t), i), cindices(t))
AbstractTrees.nodevalue(t::SpatialTree) = radixshft(first(encpoints(t)), depth(t), bitlen(t))
#!SECTION

"""
$(SIGNATURES)

Linear searches the next node high.

# Arguments
  - `R`: The array to search.
  - `lo`: The start index.
  - `hi`: The stop index.
  - `depth`: The current depth in the tree. (Starts at 0)
  - `bitlen`: The bitlen of the each bit group.
"""
function get_node_hi(R, lo, hi, depth, bitlen)
  tag(x) = radixshft(x, depth, bitlen)

  lo_tag = tag(@inbounds R[lo])
  lo += 1
  while lo <= hi && lo_tag == tag(@inbounds R[lo])
    lo += 1
  end
  return lo - 1
end


"""
$(SIGNATURES)

Counts all nodes in the tree.

# Arguments
  - `R`: The array to search.
  - `lo`: The start index.
  - `hi`: The stop index.
  - `l`: Maximum depth.
  - `depth`: The current depth in the tree. (Starts at 0)
  - `bitlen`: The bitlen of the each bit group.
"""
function count_nodes(R, lo, hi, l, depth, bitlen)
  lo > hi    && return hi, 0
  depth == l && return get_node_hi(R, lo, hi, depth, bitlen), 1

  tag(x)  = radixshft(x, depth, bitlen)

  lo_tag = tag(@inbounds R[lo])
  count = 0

  while lo <= hi && lo_tag == tag(@inbounds R[lo])
    lo, cnt = count_nodes(R, lo, hi, l, depth+1, bitlen)
    # TODO: cnt + 1 is wasteful
    count += cnt + 1
    lo += 1
  end

  return lo-1, count
end


"""
$(SIGNATURES)

Creates a tree representation of a Morton Array.

# Arguments
  - `V`: Cloud of points.
  - `R`: The array to convert.
  - `maxdpt`: Maximum depth.
  - `smlth`: Small threshold.
  - `scale`: Scalar for the original coordinates.
  - `offset`: Per dimension offset.

# Keyword Arguments
  - `dims`: Leading dimension
"""
function make_tree(V, R, I, maxdpt, smlth, bitlen, scale, offset; dims)
  _, nodes_len = count_nodes(R, 1, length(R), maxdpt, 0, bitlen)
  nodes    = Vector{NodeInfo}(undef, nodes_len)
  children = [UInt[] for _ in 1:nodes_len]

  info = TreeInfo(V, R, I, nodes, children, bitlen, maxdpt, smlth, scale, offset; dims=dims)
  @inbounds nodes[1] = NodeInfo(1, length(R), 0, 1, 0)

  tree = SpatialTree(info, 1)
  make_tree_impl(tree, 2)

  return tree
end


"""
$(SIGNATURES)

Creates a tree representation of a Morton Array.

# Arguments
  - `node`: The current node to build.
  - `idx`: Current node index.
"""
function make_tree_impl(node, idx)
  (isdeep(node) || issmall(node)) && return idx

  tag(x) = radixshft(x, depth(node)+1, bitlen(node))

  root = getroot(node)

  lo = low(node)
  while lo <= high(node)
    hi = get_node_hi(encpoints(root), lo, high(node), depth(node)+1, bitlen(node))
    @inbounds TreeInfo(root).nodes[idx] = NodeInfo(lo, hi, depth(node)+1, idx, nindex(node))
    push!(cindices(node), idx)
    lo = hi+1
    idx += 1
  end

  for child in children(node)
    idx = make_tree_impl(child, idx)
  end

  return idx
end


scalar(tree::SpatialTree) = tree.info.scale
translation(tree::SpatialTree) = tree.info.offset

qcenter!(center, node::SpatialTree) = qcenter!(center, nodevalue(node), bitlen(node), depth(node))
qcenter(T, node::SpatialTree) = qcenter!(Vector{T}(undef, bitlen(node)), node)
qcenter(node::SpatialTree) = qcenter(eltype(node), node)


function qcenter!(center, tag, bitlen, depth)
  fill!(center, 0.5)

  for _ in 1:depth, bit in 1:bitlen
    x = @inbounds center[bit]
    if tag & 1 != 0
      x += 1
    end
    @inbounds center[bit] = x/2
    tag >>>= 1
  end

  return center
end

function center!(center, node::SpatialTree)
  center .= qcenter!(center,  node) ./ scalar(node) .+ translation(node)
  return center
end
center(T, node::SpatialTree) = center!(Vector{T}(undef, bitlen(node)), node)
center(node::SpatialTree) = center(eltype(node), node)


qbox(T, node::SpatialTree)  = qbox(T, depth(node))
qbox(node::SpatialTree) = qbox(eltype(node), depth(node))

qbox(T, depth) = one(T)/(1<<(depth+1))
qbox(depth) = qbox(typeof(eps()), depth)


function box!(box, node::SpatialTree)
  box .= qbox(eltype(box), node) ./ scalar(node)
  return box
end
box(T, node::SpatialTree) = box!(Vector{T}(undef, bitlen(node)), node)
box(node::SpatialTree) = box(eltype(node), node)

Base.@propagate_inbounds function original_perm!(tree::SpatialTree, I, D)
  staticselectdim(I, Val(leaddim(tree.info)), tree.info.perm) .= tree.info.perm[I]
  staticselectdim(D, Val(leaddim(tree.info)), tree.info.perm) .= D
  return I, D
end

Base.@propagate_inbounds original_perm(tree, I, D) = original_perm!(tree, copy(I), copy(D))
