const N_REALY_BIG_ARRAY = 10000
const N_BIG_ARRAY = 100
const L_TH = Inf
const N_SMALL = 1



# NOTE: After ping-pong sorting use P array to find the true result. `Ra[P .== 1] = Rb[P .== 1]`


radixsort_seq_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, lo, hi, l, a) = @inbounds @views begin

  hi-lo+1 <= N_SMALL   && return
  l > L_TH             && return

  # TODO: Type correctness between C, lo and hi
  C = alloc!(a, 256)

  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, lo, hi, l)
  P[lo:hi] .= l%2

  for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : hi

    radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, nlo, nhi, l+1, a)
  end

  free!(a, C)
end


radixsort_par_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, lo, hi, l, a) = @inbounds @views begin

  C = alloc!(a, 256)
  countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, lo, hi, l)
  P[lo:hi] .= l%2

  @sync for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : hi

    nb = nhi-nlo+1

    nb <= 1 && continue

    if nb > N_BIG_ARRAY
      Threads.@spawn radixsort_par_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1, a)
    else
      Threads.@spawn radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1, a)
    end
  end

  free!(a, C)
end

radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, lo, hi, l, a) = @inbounds @views begin

  Cm = alloc!(a, 256, Threads.nthreads())
  C  = alloc!(a, 256)
  countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, Cm, lo, hi, l)
  minimum!(reshape(C, :, 1), Cm)
  free!(a, Cm)

  @sync for j in 1:length(C)
    nlo = C[j]+1
    nhi = j != length(C) ? C[j+1] : hi

    nb = nhi-nlo+1

    nb <= 1 && continue

    if nb > N_REALY_BIG_ARRAY
      Threads.@spawn radixsort_par_par_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1, a)
    elseif nb > N_BIG_ARRAY
      Threads.@spawn radixsort_par_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1, a)
    else
      Threads.@spawn radixsort_seq_seq_impl!(Vb, Rb, Ib, Va, Ra, Ia, P, $nlo, $nhi, l+1, a)
    end
  end

  free!(a, C)
end


