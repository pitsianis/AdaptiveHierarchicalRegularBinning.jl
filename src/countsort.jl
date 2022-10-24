using Static

countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C) = @inbounds @views begin
  fill!(C, 0)

  for r in Ra
    C[r] += 1
  end

  cumsum!(C, C)

  for i in eachindex(Ra)
    j = C[Ra[i]]
    C[Ra[i]] -= 1

    Vb[:, j] .= Va[:, i]
    Rb[j]     = Ra[i]
    Ib[j]     = Ia[i]
  end
end

countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, C) = @inbounds @views begin
  n  = length(Ra)
  np = Threads.nthreads()
  b  = cld(n, np)

  @inline local_range(k) = (k-1)*b + 1 : min(k*b, n)

  fill!(C, 0)

  Threads.@threads for k=1:np
    for i in local_range(k)
      C[Ra[i], k] += 1
    end
  end

  cumsum!(C, C, dims=2)
  cumsum!(C[:, end], C[:, end])
  C[2:end, 1:end-1] .+= C[1:end-1, end]

  Threads.@threads for k=1:np
    for i in local_range(k)
      j = C[Ra[i], k]
      C[Ra[i], k] -= 1

      Vb[:, j] .= Va[:, i]
      Rb[j]     = Ra[i]
      Ib[j]     = Ia[i]
    end
  end
end
