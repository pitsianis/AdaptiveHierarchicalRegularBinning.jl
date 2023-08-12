using AdaptiveHierarchicalRegularBinning, AbstractTrees
using BenchmarkTools

using Transducers, LoopVectorization, ThreadsX

function fastpermute_par!(O::To, I::Ti, p::Tp) where {To<:AbstractMatrix, Ti<:AbstractMatrix, Tp<:AbstractVector}

  @inbounds Threads.@threads for j in eachindex(p)
    l = p[j]
    for i in axes(O,1)
      O[i, j] = I[i, l]
    end
  end

end

function fastpermute_par!(O::To, I::Ti, p::Tp) where {To<:AbstractVector, Ti<:AbstractVector, Tp<:AbstractVector}

  @inbounds Threads.@threads for j in eachindex(p)
    l = p[j]
    O[j] = I[l]
  end

end

function countandpermute_par!(ix::Ti, Vp::Tvout, Rp::Trout, V::Tvin, R::Trin) where {Ti<:AbstractVector, Tvout<:AbstractMatrix, Trout<:AbstractVector, Tvin<:AbstractMatrix, Trin<:AbstractVector}

  sort!(ix; order = Base.Sort.Perm(Base.Order.Forward, vec(R)), alg=ThreadsX.MergeSort )
  fastpermute_par!(Rp, R, ix)
  fastpermute_par!(Vp, V, ix)

end

function fastpermute!(O::To, I::Ti, p::Tp) where {To<:AbstractMatrix, Ti<:AbstractMatrix, Tp<:AbstractVector}

  @inbounds for j in eachindex(p)
    l = p[j]
    for i in axes(O,1)
      O[i, j] = I[i, l]
    end
  end

end

function fastpermute!(O::To, I::Ti, p::Tp) where {To<:AbstractVector, Ti<:AbstractVector, Tp<:AbstractVector}

  @inbounds for j in eachindex(p)
    l = p[j]
    O[j] = I[l]
  end

end

function countandpermute!(ix::Ti, Vp::Tvout, Rp::Trout, V::Tvin, R::Trin) where {Ti<:AbstractVector, Tvout<:AbstractMatrix, Trout<:AbstractVector, Tvin<:AbstractMatrix, Trin<:AbstractVector}

  sort!(ix; order = Base.Sort.Perm(Base.Order.Forward, vec(R)) )
  fastpermute!(Rp, R, ix)
  fastpermute!(Vp, V, ix)

end

function do_interleave_parallel!(R::Vector{T}, V::Matrix{Float64}, maxL::Int, d::Int) where {T<:Unsigned}
  @inbounds @fastmath Threads.@threads for j = 1 : length(R)
    c = zero( T )
    for i = 1 : d
      v = floor( T, V[i,j] )
      for l = 1:maxL
        b = ( (T(1) << (l-1)) & v ) # find correct bit for level
        b = b << ( i-1 + (d-1)*(l-1) ) # shift to correct position
        c |= b
      end
    end
    R[j] = c
  end
end

function do_interleave_sequential!(R::Vector{T}, V::Matrix{Float64}, maxL::Int, d::Int) where {T<:Unsigned}
  @inbounds @fastmath  for j = 1 : length(R)
    c = zero( T )
    for i = 1 : d
      v = floor( T, V[i,j] )
      for l = 1:maxL
        b = ( (T(1) << (l-1)) & v ) # find correct bit for level
        b = b << ( i-1 + (d-1)*(l-1) ) # shift to correct position
        c |= b
      end
    end
    R[j] = c
  end
end

function do_interleave_parallel_old!(R::Vector{T}, V, maxL, d) where {T<:Unsigned}
  @inbounds Threads.@threads for j in axes( V, 2 )
    c = zero( T )
    for i in axes( V, 1 )
      v = floor( T, V[i,j] )
      for l = 1:maxL
        b = ( (T(1) << (l-1)) & v ) # find correct bit for level
        b = b << ( i-1 + (d-1)*(l-1) ) # shift to correct position
        @assert b & c == 0 # make sure we do not override a bit
        c |= b
      end
    end
    R[j] = c
  end
end

function do_interleave_sequential_old!(R::Vector{T}, V, maxL, d) where {T<:Unsigned}
  @inbounds for j in axes( V, 2 )
    c = zero( T )
    for i in axes( V, 1 )
      v = floor( T, V[i,j] )
      for l = 1:maxL
        b = ( (T(1) << (l-1)) & v ) # find correct bit for level
        b = b << ( i-1 + (d-1)*(l-1) ) # shift to correct position
        @assert b & c == 0 # make sure we do not override a bit
        c |= b
      end
    end
    R[j] = c
  end
end

function normalize_sequential_old!(R, V, maxL, maxdepth=maxL)
  offset = vec( minimum(V, dims=2) )
  V .-= offset
  scale  = maximum(V) + eps( eltype(V)(2)^maxdepth )
  
  V ./= scale
  V .*= eltype(R)(2)^maxL
  
  return offset, scale
end

function normalize_sequential!(R, V, maxL, maxdepth=maxL)
  ex     = foldl( TeeRF(min,max), eachcol(V) |> Broadcasting() )
  offset = map( first, ex )
  scale  = maximum( x -> x[2]-x[1], ex ) + eps( eltype(V)(2)^maxdepth )
  mult   = ( eltype(R)(2)^maxL ) / scale
  @turbo for j in axes( V, 2 )
    for i in axes(V, 1)
      V[i,j] = mult * ( V[i,j] - offset[i] )
    end
  end
  
  return offset, scale
end

function normalize_parallel!(R, V, maxL, maxdepth=maxL)
  ex     = foldl( TeeRF(min,max), eachcol(V) |> Broadcasting() )
  offset = map( first, ex )
  scale  = maximum( x -> x[2]-x[1], ex ) + eps( eltype(V)(2)^maxdepth )
  mult   = ( eltype(R)(2)^maxL ) / scale
  @tturbo for j in axes( V, 2 )
    for i in axes(V, 1)
      V[i,j] = mult * ( V[i,j] - offset[i] )
    end
  end
  
  return offset, scale
end

function compare(d, n, maxL, enctype)

  X = rand(d, n)
  
  Vgold = copy(X)
  Rgold = zeros(enctype, n)
  normalize_sequential_old!(Rgold, Vgold, maxL)

  Rigold = copy(Rgold)
  do_interleave_sequential_old!(Rigold, Vgold, maxL, d)

  Rigold_p = copy(Rigold)
  Vgold_p  = copy(Vgold)
  ixgold   = [LinearIndices(Rigold);]

  countandpermute!( ixgold, Vgold_p, Rigold_p, Vgold, Rigold )

  t_cp_seq = @benchmark( countandpermute!( ix, V_p, Ri_p, $Vgold, $Rigold ),
    setup = begin 
      if issorted($Rigold)
        error("Input should not be already sorted!")
      end
      Ri_p = copy($Rigold)
      V_p  = copy($Vgold)
      ix   = [LinearIndices($Rigold);]
    end,
    teardown = begin 
      if !( Ri_p == $Rigold_p && V_p ≈ $Vgold_p && ix == $ixgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  t_cp_par = @benchmark( countandpermute_par!( ix, V_p, Ri_p, $Vgold, $Rigold ),
    setup = begin 
      if issorted($Rigold)
        error("Input should not be already sorted!")
      end
      Ri_p = copy($Rigold)
      V_p  = copy($Vgold)
      ix   = [LinearIndices($Rigold);]
    end,
    teardown = begin 
      if !( Ri_p == $Rigold_p && V_p ≈ $Vgold_p && ix == $ixgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  t_il_seq = @benchmark( do_interleave_sequential!(R, V, $maxL, $d),
    setup = begin 
      R = copy($Rgold)
      V = copy($Vgold)
    end,
    teardown = begin 
      if !( R == $Rigold && V ≈ $Vgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  t_il_par = @benchmark( do_interleave_parallel!(R, V, $maxL, $d),
    setup = begin 
      R = copy($Rgold)
      V = copy($Vgold)
    end,
    teardown = begin 
      if !( R == $Rigold && V ≈ $Vgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )


  t_par = @benchmark( normalize_parallel!(R, V, $maxL),
    setup = begin 
      R = zeros($enctype, $n)
      V = copy($X)
    end,
    teardown = begin 
      if !( R == $Rgold && V ≈ $Vgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  t_seq = @benchmark( normalize_sequential!(R, V, $maxL),
    setup = begin 
      R = zeros($enctype, $n)
      V = copy($X)
    end,
    teardown = begin 
      if !( R == $Rgold && V ≈ $Vgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  t_old = @benchmark( normalize_sequential_old!(R, V, $maxL),
    setup = begin 
      R = zeros($enctype, $n)
      V = copy($X)
    end,
    teardown = begin 
      if !( R == $Rgold && V ≈ $Vgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  t_il_seq_old = @benchmark( do_interleave_sequential_old!(R, V, $maxL, $d),
    setup = begin 
      R = copy($Rgold)
      V = copy($Vgold)
    end,
    teardown = begin 
      if !( R == $Rigold && V ≈ $Vgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  t_il_par_old = @benchmark( do_interleave_parallel_old!(R, V, $maxL, $d),
    setup = begin 
      R = copy($Rgold)
      V = copy($Vgold)
    end,
    teardown = begin 
      if !( R == $Rigold && V ≈ $Vgold ) 
        error("Incorrect result!")
      end
    end, evals = 1 
  )

  return (;t_seq, t_par, t_old, t_il_seq_old, t_il_par_old, t_il_seq, t_il_par, t_cp_seq, t_cp_par)

end

res = compare(5, 4_000_000, 25, UInt128);

@info "Execution times (milliseconds) with $(Threads.nthreads()) threads" (map( x -> minimum(x.times) / 1e6, res ))...