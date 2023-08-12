function ahrb_fixed_length(V, maxL, maxP; QT = UInt, ctxtype = Nothing, gtctype = Nothing)::SpatialTree

  # Check that maxL is less than or equal to 25 (to avoid issues with floating-point precision)
  @assert maxL <= 25 "maxL must be <= 25"

  # Create a vector R of type QT with the same number of elements as the number of columns in V
  R = Vector{QT}(undef, size(V, 2))

  # Calculate the bit length of V
  bitlen = size(V, 1)

  # Check that the product of maxL and bitlen is less than or equal to the number of bits in QT
  maxL * bitlen > sizeof(QT) * 8 && throw("Not enough bits to represent the requested tree")

  # Create a copy of V called Vs that is going to be permuted and stored with the matrix
  Vs = copy(V)

  # Call the fast_spatial_encode! function to encode the data points in Vs into the vector R
  offset, scale = fast_spatial_encode!(R, Vs, maxL)

  # Create a linear index array ix and copy V and R into new arrays Vp and Rp
  ix = [LinearIndices(R);]
  Vp = copy(V)
  Rp = copy(R)

  # Call the countandpermute! function to bin & permute the data points in Vp and Rp
  countandpermute!(ix, Vp, Rp, V, R)

  # Call the make_tree function to generate the tree structure using 
  # the binned data points in Vp and Rp
  tree = make_tree(Vp, Rp, ix, maxL, maxP, scale, offset; ctxtype, gtctype)

  return tree

end
