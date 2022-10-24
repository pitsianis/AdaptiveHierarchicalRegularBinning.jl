using AdaptiveHierarchicalRegularBinning: countsort_seq_impl!, countsort_par_impl!
using Test

@testset "countsort" begin
  n = 10_000
  d = 3

  Va = rand(d, n)
  Ra = rand(UInt8, n)
  Ia = UInt.(1:n)

  Vb = similar(Va)
  Rb = similar(Ra)
  Ib = similar(Ia)

  @testset "Sequential" begin
    C = Vector{UInt}(undef, maximum(Ra)+1)
    countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C)

    @test issorted(Rb)
    @test Ra == Rb[sortperm(Ib)]
    @test Va == Vb[:, sortperm(Ib)]
  end

  @testset "Parallel" begin
    C = Matrix{UInt}(undef, maximum(Ra)+1, Threads.nthreads())
    countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, C)

    @test issorted(Rb)
    @test Ra == Rb[sortperm(Ib)]
    @test Va == Vb[:, sortperm(Ib)]
  end
end
