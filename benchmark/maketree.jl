using AdaptiveHierarchicalRegularBinning
using BenchmarkTools

function doit(n,d,l)
 
  X = rand(d, n)

  @btime tree = regular_bin(UInt128, $X, $l, 2^8; dims=2)

  nothing
end
