using AdaptiveHierarchicalRegularBinning: countsort_seq_impl!, countsort_par_impl!, CountSortDetails, bitmask, radixsel
using Test


n = 10_000
d = 3
@testset "CountSort(n=$n, d=$d)" begin

  @testset "CountSortDetails(bitlen=$bitlen, leaddim=$dims)" for bitlen = 1:16, dims = 1:2
    print("bitlen=$bitlen, leaddim=$dims")

    Va = dims == 1 ? rand(n, d) : rand(d, n)
    Ra = rand(UInt, n)
    Ia = UInt.(1:n)

    Vb = similar(Va)
    Rb = similar(Ra)
    Ib = similar(Ia)


    @testset "Sequential" begin
      csd = CountSortDetails(bitlen, 1, n; dims=dims)
      C = Vector{UInt}(undef, bitmask(csd)+1)
      countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, csd)

      @test issorted( radixsel.(Ref(csd),  Rb) )
      @test Ra == Rb[sortperm(Ib)]
      @test Va == (dims == 1 ? Vb[sortperm(Ib), :] : Vb[:, sortperm(Ib)])

      lo = UInt(1 + n÷4)
      hi = UInt(n÷2)
      csd = CountSortDetails(bitlen, lo, hi; dims=dims)
      Rb .= Ra
      countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, csd)

      @test issorted( radixsel.(Ref(csd),  Rb[lo:hi]) )
      print(" Sequential OK.")
    end


    @testset "Parallel" begin
      csd = CountSortDetails(bitlen, 1, n; dims=dims)
      C = Matrix{UInt}(undef, bitmask(csd)+1, Threads.nthreads())
      countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, csd)

      @test issorted( radixsel.(Ref(csd),  Rb) )
      @test Ra == Rb[sortperm(Ib)]
      @test Va == (dims == 1 ? Vb[sortperm(Ib), :] : Vb[:, sortperm(Ib)])

      lo = UInt(1 + n÷4)
      hi = UInt(n÷2)
      csd = CountSortDetails(bitlen, lo, hi; dims=dims)
      Rb .= Ra
      countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, csd)

      @test issorted( radixsel.(Ref(csd),  Rb[lo:hi]) )
      println(" Parallel OK.")
    end
  end

end
