using AdaptiveHierarchicalRegularBinning: countsort_seq_impl!, countsort_par_impl!
using BenchmarkTools
using Plots



function bmark_single(countsort_impl!, d, n, TR=UInt8)
  B = @benchmark(
    $countsort_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, lo, hi, nbyte),
    setup=begin
      Va = rand($d, $n)
      Ra = rand($TR, $n)
      Ia = UInt.(1:$n)

      Vb = similar(Va)
      Rb = similar(Ra)
      Ib = similar(Ia)

      lo = UInt(1)
      hi = UInt(length(Ra))

      nbyte = UInt(1)

      C = Dict(
        countsort_seq_impl! => ()->Vector{UInt}(undef, maximum(Ra)+1),
        countsort_par_impl! => ()->Matrix{UInt}(undef, maximum(Ra)+1, Threads.nthreads())
      )[$countsort_impl!]()
    end
  )

  return median(B.times)
end


d = 10
N = 10 .^ (3:6)

R = Dict()

R["seq"] = []
for n in N
  append!(R["seq"], bmark_single(countsort_seq_impl!, d, n))
end

R["par"] = []
for n in N
  append!(R["par"], bmark_single(countsort_par_impl!, d, n))
end


plot()
plot!(N, R["seq"]; label="Sequential")
plot!(N, R["par"]; label="Parallel")
