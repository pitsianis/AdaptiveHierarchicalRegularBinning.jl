using AdaptiveHierarchicalRegularBinning: radixsort_seq_seq_impl!, radixsort_par_seq_impl!, radixsort_par_par_impl!, RadixSortDetails, Allocator
using BenchmarkTools

n = 2_000_000
d = 8
l = 8
sigma = 2.0

@inbounds function do_benchmark(n, d, l, sigma)
  V = sigma * randn(d, n)
  R = zeros(UInt, n)
  AdaptiveHierarchicalRegularBinning.spatial_encode!(R,V,l; dims = Val(2), center = false)
  I = collect(UInt, 1:n)

  # ("Seq2Seq", "Par2Seq", "Par2Par"), 
  # (radixsort_seq_seq_impl!, radixsort_par_seq_impl!, radixsort_par_par_impl!))

  Va = copy(V);    Ra = copy(R);    Ia = copy(I)
  Vb = similar(V); Rb = similar(R); Ib = similar(I)

  P = zeros(Bool, n)

  rsd = RadixSortDetails(d, 1, n; dims=2, sml_th=1, dpt_th=l + 1)

  allocator = Allocator(eltype(I))

  println("Par2Par")
  @btime begin

    radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, allocator)
    Va[:, P] .= Vb[:, P]
    Ra[P] .= Rb[P]
    Ia[P] .= Ib[P]
  
  end setup = (Va = deepcopy($Va); Ra = deepcopy($Ra); Ia = deepcopy($Ia);
               Vb = deepcopy($Vb); Rb = deepcopy($Rb); Ib = deepcopy($Ib);
               P = deepcopy($P); rsd = deepcopy($rsd); allocator = deepcopy($allocator))

  radixsort_par_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, allocator)
  Va[:, P] .= Vb[:, P]
  Ra[P] .= Rb[P]
  Ia[P] .= Ib[P]

  pVa = similar(V)
  println("sortperm")
  @btime begin
    pIa = sortperm(R)
    pVa .= $V[:, pIa]
  end setup=( R = copy($R); V = copy($V); pVa = copy($pVa) )

  pIa = sortperm(R)
  pVa .= V[:, pIa]

  @assert pIa == Ia
  @assert pVa == Va

end

do_benchmark(n, d, l, sigma)