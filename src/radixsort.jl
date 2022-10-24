const N_REALY_BIG_ARRAY = 10000
const N_BIG_ARRAY = 100
const L_TH = Inf
const N_SMALL = 1



# NOTE: After ping-pong sorting use P array to find the true result. `Ra[P .== 1] = Rb[P .== 1]`


radixsort_seq_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, lo, hi, l) = @inbounds @views begin

  hi-lo+1 <= N_SMALL   && return
  l > L_TH             && return

  # TODO: Allocator
  C = Vector{UInt}(undef, 256)

  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, lo, hi, l)
  P[lo:hi] .= l%2

  for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : hi

    radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, nlo, nhi, l+1)
  end
end


radixsort_par_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, lo, hi, l) = @inbounds @views begin
  C = Vector{UInt}(undef, 256)
  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, lo, hi, l)
  P[lo:hi] .= l%2

  @sync for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : hi

    nb = nhi-nlo+1

    # TODO: Early continue for nb <= 1

    if nb > N_BIG_ARRAY
      Threads.@spawn radixsort_par_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1)
    else
      Threads.@spawn radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1)
    end
  end

end

radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, lo, hi, l) = begin

  Cm = Matrix{UInt}(undef, 256, Threads.nthreads())
  C  = Cm[:, end]
  countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, Cm, lo, hi, l)
  minimum!(C, Cm)

  for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : hi

    nb = nhi-nlo+1

    if nb > N_REALY_BIG_ARRAY
      Threads.@spawn radixsort_par_par_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1)
    elseif nb > N_BIG_ARRAY
      Threads.@spawn radixsort_par_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1)
    else
      Threads.@spawn radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1)
    end
  end
end


