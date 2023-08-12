
"""
    $(TYPEDSIGNATURES)

Find leaves that are above the maximum number of points per leaf and we should recurse to subdivide
them.
"""
function find_rec_leaves(nodes::Vector{NodeInfo}, children::Vector{Vector{Int64}}, maxpts::Int64)

  filter( i -> isempty(children[i]) && (nodes[i].hi - nodes[i].lo + 1) > maxpts, 1:length(nodes) );

end

"""
    $(TYPEDSIGNATURES)

Copy the codes and the permutation vector between the ping-pong buffers for the final result
"""
function copy_final!(R0::Tr, R1::Tr, data_0or1::Vector{Int8}) where {Tr<:AbstractVector}
  
  @assert length(R0) == length(R1) == length(data_0or1)
  
  @inbounds for i = eachindex(R0)
    if data_0or1[i] == 1
      R0[i]  = R1[i]
    end
  end

end

"""
    $(TYPEDSIGNATURES)

Compute the center and sidelength of each node in the tree.
"""
function compute_centers_sidelengths!(centers::Vector{Vector{V}}, sidelengths::Vector{V}, 
    nodes::Vector{NodeInfo}, R0::Vector{Tu}, currentdepth::Int, d::Int) where {V<:Real, Tu<:Unsigned}
  
  @inbounds for i = eachindex( nodes )
    lo  = nodes[i].lo
    dpt = nodes[i].dpt
    _qcenter!(centers[i], R0[lo], d, dpt, currentdepth)
    sidelengths[i] = sidelengths[i] * 2.0^(-dpt)
  end

end

"""
    $(TYPEDSIGNATURES)

Select which buffer to use for the ping-pong algorithm.
"""
function select_pingpong_buffers_subdivision(R0, V0, Pglb, R1, V1, Plcl, rng, data_0or1)
  
  @inbounds begin
    if data_0or1 == 0  # data in 0, we write to 1
      Rout = @view R1[rng]
      Vout = @view V1[:, rng]

      Rin  = @view R0[rng]
      Vin  = @view V0[:, rng]
    else
      Rout = @view R0[rng]
      Vout = @view V0[:, rng]

      Rin  = @view R1[rng]
      Vin  = @view V1[:, rng]
    end
  end

  Plcl_sub = @view Plcl[ rng ]      
  Pglb_sub = @view Pglb[ rng ]

  return Rout, Vout, Plcl_sub, Rin, Vin, Pglb_sub

end

@inbounds function ahrb_block(V, maxdepth, maxpoints; lstep = 2, QT = UInt, ctxtype = Nothing, gtctype = Nothing)

  @assert maxdepth <= 25 "deepestlevel must be <= 25"

  # Total number of particles
  (d, n) = size(V)

  # Create a vector R of type QT with the same number of elements as the number of columns in V
  R0 = zeros(QT, size(V, 2))

  # Calculate the bit length of V
  bitlen = size(V, 1)

  # Check that the product of lstep and bitlen is less than or equal to the number of bits in QT
  lstep * bitlen > sizeof(QT) * 8 && throw("Not enough bits to represent the requested tree")

  ## --- First lstep levels processing --- ##

  # Create a copy of V called Vs that is going to be permuted and stored with the matrix
  V0 = copy(V)

  # Call the fast_spatial_encode! function to encode the data points in Vs into the vector R
  offset, scale = fast_spatial_encode!(R0, V0, lstep, maxdepth)

  # Create a linear index array ix and copy V and R into new arrays Vp and Rp
  Pglb = [LinearIndices(R0);]
  Plcl = copy(Pglb)
  V1   = copy(V0)
  R1   = copy(R0)

  # Call the countandpermute! function to bin & permute the data points in Vp and Rp
  countandpermute!(Pglb, V1, R1, V0, R0)

  # The current maximum level is the lstep
  currentdepth = lstep

  # Build the nodes and the children pointers only for the first lstep levels
  nodes = NodeInfo[]
  children = Vector{Vector{Int64}}()
  sizehint!(nodes, 2*n)
  sizehint!(children, 2*n)

  # Add the root to the tree
  push!( nodes, NodeInfo(1, length(R1), 0, 1, 0) )
  push!( children, Int64[] )
  
  # Form the tree up to this level
  form_tree!( nodes, children, nodes[1], R1, currentdepth, maxpoints, d )
  update_children!( nodes, children )

  ## --- Subsequent levels processing recursively --- ##

  # Now recursively subdivide leaves until maxP is reached or 25 levels are reached, or the bits are
  # exhausted
  idx = find_rec_leaves( nodes, children, maxpoints )
  
  data_0or1 = ones(Int8, n)  # 0 or 1 buffer

  while !isempty(idx) && ( currentdepth < maxdepth ) && ( (currentdepth + lstep) * bitlen <= sizeof(QT) * 8 )

    # Shift the bits of all particles for new level (prepare space for new bits)
    prepare_next_levels!(R0, R1, data_0or1, lstep, d)
    currentdepth += lstep

    # Subdivide leaves that are above maxP
    @threading for i in idx
      partidx = range( nodes[i] )
      which_buffer = data_0or1[ partidx[1] ]

      # Select the correct parts of the ping-pong buffers for this block
      Rout, Vout, Plcl_sub, Rin, Vin, Pglb_sub = select_pingpong_buffers_subdivision(
        R0, V0, Pglb, R1, V1, Plcl, partidx, which_buffer )
      
      # initialize local permutation vector
      Plcl_sub .= 1:length(partidx)

      # local encoding of the new bits
      fast_spatial_encode_block!( Rin, Vin, lstep )

      # local binning and permutation
      countandpermute_seq!(Plcl_sub, Vout, Rout, Vin, Rin)

      @debugging begin
        @assert Vin[ :, Plcl_sub ] == Vout
        @assert Rin[ Plcl_sub ] == Rout
        @assert isperm( Plcl_sub )
        @assert minimum( Plcl_sub ) == 1 && maximum( Plcl_sub ) == length(partidx)
      end

      # update on which buffer the data for this node & its children is stored
      data_0or1[partidx] .= 1 - which_buffer

      # update the global permutation vector
      Pglb_sub .= Pglb_sub[Plcl_sub]

    end

    for i in idx
      # find where to build this tree
      which_buffer = data_0or1[ nodes[i].lo ]
      Rout = which_buffer == 0 ? R0 : R1
      
      # form the subtree at this node
      idx_new = form_tree!( nodes, children, nodes[i], Rout, currentdepth, maxpoints, d )
      update_children!( @view( nodes[idx_new] ), children )
    end

    # Find the leaves that need to be subdivided for the next iteration
    idx = find_rec_leaves( nodes, children, maxpoints )

  end

  ## --- Final processing --- ##

  # copy the data from 1 to 0 for final processing
  copy_final!(R0, R1, data_0or1)

  # Compute center and sidelength for each box
  centers = [zeros(eltype(V),d) for _ = 1 : length(nodes)]
  sidelengths = ones(eltype(V), length(nodes))
  compute_centers_sidelengths!(centers, sidelengths, nodes, R0, currentdepth, d)

  # Build tree data
  info = TreeInfo(V[:,Pglb], R0, Pglb, nodes, children, currentdepth, maxpoints, scale, offset, centers, sidelengths; ctxtype, gtctype)
  resize!( info.context, length(nodes) )

  tree = SpatialTree{eltype(V),QT,ctxtype,gtctype}(info, 1)

  return tree

end