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
  - `scale`: The original scale of the dimensions of the `points`.
  - `offset`: The per dimension offset of the `points`.
  - `nodes`: Per node information.
  - `children`: Per node children.
  - `maxdepth`: The maximum depth of the tree.
"""
struct TreeInfo{T, E, B, D}
  points::Matrix{T}
  encoded::Vector{E}

  scale::T
  offset::Vector{T}

  nodes::Vector{NodeInfo}
  children::Vector{Vector{UInt}}

  maxdepth::UInt
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
TreeInfo(V, R, N, C, bitlen, depth, scale, offset; dims) = TreeInfo{eltype(V), eltype(R), bitlen, dims}(V, R, scale, offset, N, C, depth)
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
  - `R`: The array to convert.
  - `l`: Maximum depth.
  - `bitlen`: The bitlen of the each bit group.
"""
function make_tree(V, R, maxdpt, bitlen, scale, offset; dims)
  _, nodes_len = count_nodes(R, 1, length(R), maxdpt, 0, bitlen)
  nodes    = Vector{NodeInfo}(undef, nodes_len)
  children = [UInt[] for _ in 1:nodes_len]

  info = TreeInfo(V, R, nodes, children, bitlen, maxdpt, scale, offset; dims=dims)
  @inbounds nodes[1] = NodeInfo(1, length(R), 0, 1, 0)

  tree = SpatialTree(info, 1)
  make_tree_impl(tree, 2)

  return tree
end


"""
$(SIGNATURES)

Creates a tree representation of a Morton Array.

# Arguments
  - `tree`: The tree information.
  - `root`: The current tree node.
  - `l`: Maximum depth.
  - `idx`: Current node index.
"""
function make_tree_impl(node, idx)
  isdeep(node) && return idx

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


qcenter!(center, node::SpatialTree) = qcenter!(center, nodevalue(node), bitlen(node), depth(node))
qcenter(T, node::SpatialTree) = qcenter(T, nodevalue(node), bitlen(node), depth(node))
qcenter(node::SpatialTree) = qcenter(eltype(node), nodevalue(node), bitlen(node), depth(node))


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
end

function qcenter(T, tag, bitlen, depth)
  center = Vector{T}(undef, bitlen)
  qcenter!(center, tag, bitlen, depth)
  return center
end

qcenter(tag, bitlen, depth) = qcenter(typeof(eps()), tag, bitlen, depth)



qbox(T, node::SpatialTree)  = qbox(T, depth(node))
qbox(node::SpatialTree) = qbox(eltype(node), depth(node))

qbox(T, depth) = one(T)/(1<<depth)
qbox(depth) = qbox(typeof(eps()), depth)

