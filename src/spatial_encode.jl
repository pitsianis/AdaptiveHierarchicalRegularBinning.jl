function do_interleave!(R::Vector{T}, V::Matrix{Float64}, maxL::Int, d::Int) where {T<:Unsigned}
  @inbounds @fastmath @threading for j = 1 : length(R)
    c = zero( T )
    for i = 1 : d
      v = floor( T, V[i,j] )
      for l = 1:maxL
        b = ( (T(1) << (l-1)) & v ) # find correct bit for level
        b = b << ( i-1 + (d-1)*(l-1) ) # shift to correct position
        # @assert b & c == 0 # make sure we do not override a bit
        c |= b
      end
    end
    R[j] = c
  end
end

function fast_spatial_encode!(R, V, maxL, maxdepth=maxL)
  ex     = foldl( TeeRF(min,max), eachcol(V) |> Broadcasting() )
  offset = map( first, ex )
  scale  = maximum( x -> x[2]-x[1], ex ) + eps( eltype(V)(2)^maxdepth )
  mult   = ( eltype(V)(2)^maxL ) / scale
  d      = size(V,1)
  @tturbo for j in axes( V, 2 )
    for i in axes(V, 1)
      V[i,j] = mult * ( V[i,j] - offset[i] )
    end
  end
  
  do_interleave!(R, V, maxL, d)

  return offset, scale
end

function do_interleave_block!(R::Tr, V::Tv, lnew::Int, d::Int) where {T<:Unsigned, Tr<:AbstractVector{T}, Tv<:AbstractMatrix{Float64}}
  @inbounds @fastmath for j in axes( V, 2 )
    c = R[j]
    for i in axes( V, 1 )
      v = floor( T, V[i,j] )
      for l = 1:lnew
        b = ( (T(1) << (l-1)) & v ) # find correct bit for level
        b = b << ( i-1 + (d-1)*(l-1) ) # shift to correct position
        # @assert b & c == 0 # make sure we do not override a bit
        c |= b
      end
    end
    R[j] = c
  end
end

function fast_spatial_encode_block!(R, V, lnew)
  d    = size(V,1)
  mult = eltype(V)(2)^lnew
  @tturbo for j in axes( V, 2 )
    for i in axes(V, 1)
      V[i,j] = mult * V[i,j]
    end
  end

  do_interleave_block!(R, V, lnew, d)

  return nothing
end

function prepare_next_levels!(R, lnew, d)

  @inbounds @threading for j in eachindex( R )
      R[j] <<= lnew * d
  end

  return nothing

end

function prepare_next_levels!(R0, R1, data_0or1, lnew, d)

  @inbounds @threading for j in eachindex( R0 )
    if data_0or1[j] == 0
      R0[j] <<= lnew * d
    elseif data_0or1[j] == 1
      R1[j] <<= lnew * d
    else
      error("data_0or1[$j] == $(data_0or1[j])")
    end
  end

  return nothing

end