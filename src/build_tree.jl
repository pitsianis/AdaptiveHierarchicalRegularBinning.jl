function make_tree(V, R, I, maxdpt, maxpoints, scale, offset; ctxtype = Bool, gtctype = Nothing)

  n = size(V, 2)
  nodes = NodeInfo[]
  children = Vector{Vector{Int64}}()
  centers = Vector{Vector{eltype(V)}}()
  sidelengths = Vector{eltype(V)}()

  sizehint!(nodes, 2*n)
  sizehint!(children, 2*n)

  info = TreeInfo(V, R, I, nodes, children, maxdpt, maxpoints, scale, offset, centers, sidelengths; ctxtype, gtctype)
  push!( nodes, NodeInfo(1, length(R), 0, 1, 0) )
  push!( children, Int64[] )
  tree = SpatialTree(info, 1)
  form_tree!(tree)
  resize!( info.context, length(nodes) )

  for i = eachindex( nodes )
    node = SpatialTree(info, i)
    push!( centers, _qcenter(node) )
    push!( sidelengths, _qsidelength(node) )
  end

  return tree

end

function decodeboxids!(boxids::T, r::UInt128, l::Int, d::Int) where {T<:AbstractVector{UInt128}}
  mask = UInt128( 2^d-1 )
  for i = l:-1:1
    boxids[i] = r & mask
    r >>= d
  end
end

@inbounds function form_tree!(node)

  nodes = node.info.nodes
  children = node.info.children

  maxlevel = Int( node.info.maxdepth )
  d = size( node.info.points, 1 )

  Rbit = node.info.encoded

  curlevel = depth(node) + 1

  idxbeg = low(node); idxfin = low(node);
  shift = (maxlevel-curlevel)*d
  
  beg = Rbit[idxbeg] >> shift; 

  for i = low(node):high(node)

    fin = i < high(node) ? Rbit[i+1] >> shift : beg+1

    if fin != beg
      idxfin = i

      newnode = NodeInfo(idxbeg, idxfin, curlevel, length(nodes)+1, nindex(node))
      push!( nodes, newnode )
      push!( children, Int64[] )

      if curlevel < maxlevel && idxfin-idxbeg+1 > node.info.maxpoints
        form_tree!( SpatialTree( TreeInfo(node), nindex(newnode) ) )
      end

      beg    = fin
      idxbeg = i+1
    end
    
  end

  if depth(node) == 0
    for n in nodes
      if n.pidx > 0
        push!( children[n.pidx], n.nidx )
      end
    end
  end

  return nothing

end

@inline function update_children!(node::NodeInfo, children::Vector{Vector{Int64}})
  @inbounds if node.pidx > 0
    push!( children[node.pidx], node.nidx )
  end
end

@inline function update_children!(nodes::Tn, children::Vector{Vector{Int64}}) where {Tn<:AbstractVector{NodeInfo}}
  @inbounds for n in nodes
    update_children!(n, children)
  end
end

@inbounds function form_tree!(nodes::Vector{NodeInfo}, children::Vector{Vector{Int64}},
    node::NodeInfo, Rbit::Vector{T}, maxlevel::Int, maxpts::Int, d::Int ) where {T<:Unsigned}

  curlevel = depth(node) + 1

  idx_first_new = length(nodes)+1

  idxbeg = low(node); idxfin = low(node);
  shift = (maxlevel-curlevel)*d
  
  beg = Rbit[idxbeg] >> shift; 

  for i = low(node):high(node)

    fin = i < high(node) ? Rbit[i+1] >> shift : beg+1

    if fin != beg
      idxfin = i

      newnode = NodeInfo(idxbeg, idxfin, curlevel, length(nodes)+1, nindex(node))
      push!( nodes, newnode )
      push!( children, Int64[] )

      if curlevel < maxlevel && idxfin-idxbeg+1 > maxpts
        form_tree!( nodes, children, newnode, Rbit, maxlevel, maxpts, d )
      end

      beg    = fin
      idxbeg = i+1
    end
    
  end

  idx_last_new = length(nodes)

  return idx_first_new:idx_last_new

end

