function fastpermute!(O::To, I::Ti, p::Tp) where {To<:AbstractMatrix, Ti<:AbstractMatrix, Tp<:AbstractVector}

  @inbounds @threading for j in eachindex(p)
    l = p[j]
    for i in axes(O,1)
      O[i, j] = I[i, l]
    end
  end

end

function fastpermute!(O::To, I::Ti, p::Tp) where {To<:AbstractVector, Ti<:AbstractVector, Tp<:AbstractVector}

  @inbounds @threading for j in eachindex(p)
    l = p[j]
    O[j] = I[l]
  end

end

function countandpermute!(ix::Ti, Vp::Tvout, Rp::Trout, V::Tvin, R::Trin) where {Ti<:AbstractVector, Tvout<:AbstractMatrix, Trout<:AbstractVector, Tvin<:AbstractMatrix, Trin<:AbstractVector}

  sort!(ix; order = Base.Sort.Perm(Base.Order.Forward, vec(R)), alg=ThreadsX.MergeSort )
  fastpermute!(Rp, R, ix)
  fastpermute!(Vp, V, ix)

end

function fastpermute_seq!(O::To, I::Ti, p::Tp) where {To<:AbstractMatrix, Ti<:AbstractMatrix, Tp<:AbstractVector}

  @inbounds for j in eachindex(p)
    l = p[j]
    for i in axes(O,1)
      O[i, j] = I[i, l]
    end
  end

end

function fastpermute_seq!(O::To, I::Ti, p::Tp) where {To<:AbstractVector, Ti<:AbstractVector, Tp<:AbstractVector}

  @inbounds for j in eachindex(p)
    l = p[j]
    O[j] = I[l]
  end

end

function countandpermute_seq!(ix::Ti, Vp::Tvout, Rp::Trout, V::Tvin, R::Trin) where {Ti<:AbstractVector, Tvout<:AbstractMatrix, Trout<:AbstractVector, Tvin<:AbstractMatrix, Trin<:AbstractVector}

  sort!(ix; order = Base.Sort.Perm(Base.Order.Forward, vec(R)) )
  fastpermute_seq!(Rp, R, ix)
  fastpermute_seq!(Vp, V, ix)

end