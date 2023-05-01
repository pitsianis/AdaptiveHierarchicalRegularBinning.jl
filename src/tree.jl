using AbstractTrees

#SECTION: Structs
#TODO: Generic types

"""
$(TYPEDEF)

Represents a node of the tree.

# Fields
  - `lo`: The start index of the represented group.
  - `hi`: The stop index of the represented group.
  - `depth`: The depth of the node
  - `idx`: The index of the node.
"""
struct Node
  lo::UInt
  hi::UInt
  depth::UInt
  idx::UInt
end


"""
$(TYPEDEF)

Represents the tree information

# Fields
  - `data`: The data vector of the tree.
  - `nodes`: The vector of nodes of the tree.
  - `children`: The children of each node.
  - `bitlen`: The bit length of the bit groups.
"""
struct TreeInfo
  data::Vector{UInt}
  nodes::Vector{Node}
  children::Vector{Vector{UInt}}
  bitlen::UInt
end


"""
$(TYPEDEF)

Represents the tree.

# Fields
  - `info`: The tree information.
  - `root`: Index to the root of the tree.
"""
struct TreeHandle
  info::TreeInfo
  root::UInt
end
#!SECTION

#SECTION: AbstractTrees
AbstractTrees.children(t::TreeHandle)  = TreeHandle.(Ref(t.info), t.info.children[t.root])
AbstractTrees.nodevalue(t::TreeHandle) = t.info.data[t.info.nodes[t.root].lo]
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
function make_tree(R, l, bitlen)
  _, nodes_len = count_nodes(R, 1, length(R), l, 0, bitlen)
  nodes     = Vector{Node}(undef, nodes_len)
  children  = Vector{Vector{UInt}}(undef, nodes_len)
  for i in 1:nodes_len
    @inbounds children[i] = []
  end

  info = TreeInfo(R, nodes, children, bitlen)
  root = Node(1, length(R), 0, 1)
  @inbounds info.nodes[1] = root
  make_tree_impl(info, root, l, 2)

  return TreeHandle(info, 1)
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
function make_tree_impl(tree, root, l, idx)
  root.depth == l && return idx

  tag(x) = radixshft(x, root.depth+1, tree.bitlen)

  lo = root.lo
  while lo <= root.hi
    hi = get_node_hi(tree.data, lo, root.hi, root.depth+1, tree.bitlen)
    tree.nodes[idx] = Node(lo, hi, root.depth+1, idx)
    push!(@inbounds(tree.children[root.idx]), idx)
    lo = hi+1
    idx += 1
  end

  for child_idx in @inbounds(tree.children[root.idx])
    idx = make_tree_impl(tree, tree.nodes[child_idx], l, idx)
  end

  return idx
end
