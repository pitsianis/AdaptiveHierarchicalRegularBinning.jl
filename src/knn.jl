function spcolmuldist!(xb, i, A, B)
  rowvalA = rowvals(A)
  nzvalA = nonzeros(A)
  rowvalB = rowvals(B)
  nzvalB = nonzeros(B)
  @inbounds begin
    for jp in nzrange(B, i)
      nzB = nzvalB[jp]
      j = rowvalB[jp]
      for kp in nzrange(A, j)
        nzC = nzvalA[kp] * nzB
        k = rowvalA[kp]
        xb[k] -= 2nzC
      end
    end
  end
  return nothing
end

function norm2_and_kselect_col!( idx_i, dst_i, xb, i, Ct, Q, q_nrmsq_i, C_idx )
  
  xb .= xb .+ q_nrmsq_i
  
  spcolmuldist!(xb, i, Ct, Q)
  
  kselect!( idx_i, dst_i, xb, C_idx )

end

function onthefly_block_knn!(Idxs, Dists, xb, Ct, Q, C_nrmsq, Q_nrmsq, C_idx)

  nC, nD = size(Ct)
  nQ = size(Q, 2)
  nD == size(Q, 1) || throw(DimensionMismatch())
  nC == length(xb) || throw(DimensionMismatch())

  maxdist = 0.0
  
  @inbounds for i in 1:nQ
    dst_i = @view Dists[:,i]
    idx_i = @view Idxs[:,i]
    xb .= C_nrmsq
    norm2_and_kselect_col!(idx_i, dst_i, xb, i, Ct, Q, Q_nrmsq[i], C_idx)
    maxdist = max(maxdist, dst_i[1])
  end

  return maxdist

end


@inbounds function kselect!(Idxs::Tm,Dists::Tm2,v::Tv,idx::Tn) where {Tm<:AbstractVector,Tm2<:AbstractVector,Tv<:AbstractVector,Tn<:AbstractVector}
  for i in 1:length(v)
    if Dists[1] >= v[i]
      Dists[1] = v[i]
      Idxs[1] = idx[i]
      percolate_down!(Dists, Idxs, v[i], idx[i])
    end
  end
end

@inline getleft(i::Int) = 2i
@inline getright(i::Int) = 2i + 1
@inline getparent(i::Int) = div(i, 2)

# In place heap sort
@inline function heap_sort_inplace!(xs, xis)
  @inbounds for i in length(xs):-1:2
      xs[i], xs[1] = xs[1], xs[i]
      xis[i], xis[1] = xis[1], xis[i]
      percolate_down!(xs, xis, xs[1], xis[1], i - 1)
  end
  return
end

# Binary max-heap percolate down.
@inline function percolate_down!(xs::AbstractArray,
                       xis::AbstractArray,
                       dist::Number,
                       index::Int,
                       len::Int=length(xs))
  i = 1
  @inbounds while (l = getleft(i)) <= len
      r = getright(i)
      j = ifelse(r > len || (xs[l] > xs[r]), l, r)
      if xs[j] > dist
          xs[i] = xs[j]
          xis[i] = xis[j]
          i = j
      else
          break
      end
  end
  xs[i] = dist
  xis[i] = index
  return
end