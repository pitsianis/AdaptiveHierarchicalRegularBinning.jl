
"""
$(TYPEDSIGNATURES)

Multithreaded count sort permutation. Overwrite vector `p` with the
permutation that puts vector of values `v[p]` in non decreasing order.

 - `p`: preallocated output vector of integers to receive the permutation
 - `v`: input vector of integer values that determine order. The unique values
        of `v` must be in the range `0:m` where `m << length(v)`
 - `C`: preallocated workspace of size at least `(m+1+8:np+1)`,
        where `np` is the number of threads.
"""
function
pcountsortperm!(p::Vector{Int64},v::Vector{Int64},C::Matrix{Int64})

  n  = length(v)
  np = Threads.nthreads()

  fill!( C, 0 )

  b = Int64( ceil( n / np ) )

  @inbounds @threading for k = 1:np
    for i=(k - 1) * b + 1 : min(k*b,n)
      x = v[i]
      C[x, k+1] += 1
    end
  end

  cumsum!(C, C, dims=2)
  C[2:end,1:end-1] .= C[2:end,1:end-1] .+ cumsum(C[1:end-1,end])

  @inbounds @threading for k=1:np
    for i=(k - 1) * b + 1 : min(k*b,n)
      C[v[i],k] += 1
      p[C[v[i],k]] = i
    end
  end
end

"""
$(TYPEDSIGNATURES)

Sequential count sort permutation. Overwrite vector `p` with the
permutation that puts vector of values `v[p]` in non decreasing order.

 - `p`: preallocated output vector of integers to receive the permutation
 - `v`: input vector of integer values that determine order. The unique values
        of `v` must be in the range `0:m` where `m << length(v)`
 - `C`: preallocated workspace of size at least `m+1`
"""
function countsortperm!(p,v,C)

  fill!( C, 0 )

  n = length(v)

  Coff = view( C, 2:length(C) )
  @inbounds for i=1:n
    Coff[v[i]] += 1
  end

  cumsum!(C, C)

  @inbounds for i=1:n
    C[v[i]] += 1
    p[C[v[i]]] = i
  end
end


countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C) = @inbounds @views begin
  fill!(C, 0)

  Coff = C[2:end]
  for r in Ra
    Coff[r] += 1
  end

  cumsum!(C, C)

  for i in eachindex(Ra)
    C[Ra[i]] += 1
    j = C[Ra[i]]

    Vb[:, j] = Va[:, i]
    Rb[j]    = Ra[i]
    Ib[j]    = Ia[i]
  end
end
