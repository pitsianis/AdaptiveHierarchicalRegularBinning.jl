using AdaptiveHierarchicalRegularBinning: radixsort_seq_seq_impl!, radixsort_par_seq_impl!, radixsort_par_par_impl!, RadixSortDetails, Allocator
using Test


@testset "RadixSort" begin
  n = 1_000_000
  d = 8
  l = 8

  V = rand(d, n)
  R = rand(UInt, n)
  I = collect(UInt, 1:n)


  @testset "$testname" for (testname, func!) in zip(("Seq2Seq", "Par2Seq", "Seq2Seq"), (radixsort_seq_seq_impl!, radixsort_par_seq_impl!, radixsort_par_par_impl!))
    Va = copy(V)
    Ra = copy(R)
    Ia = copy(I)

    Vb = similar(V)
    Rb = similar(R)
    Ib = similar(I)

    P = falses(n)

    rsd = RadixSortDetails(d, 1, n; dims=2, sml_th=1, dpt_th=l+1)

    allocator = Allocator(eltype(I))

    func!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, allocator)

    Va[:, P] .= Vb[:, P]
    Ra[P]    .= Rb[P]
    Ia[P]    .= Ib[P]

    @test issorted(Ra)
    @test isperm(Ia)
    @test R == Ra[sortperm(Ia)]
    @test V == Va[:, sortperm(Ia)]
  end
end
