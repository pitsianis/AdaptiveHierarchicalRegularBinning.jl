using AdaptiveHierarchicalRegularBinning: countsortperm!, pcountsortperm!
using Test

@testset "seq and parallel count sort" begin
  n = 40_000
  m = 7
  np = Threads.nthreads()
  v = Int64.(ceil.(rand(n) .* m))
  p = similar(v)
  q = similar(v)
  Caux = zeros(Int64,m+8,np+1)
  C    = zeros(Int64,m+1)

  pcountsortperm!(p,v,Caux)
  countsortperm!(q,v,C)

  @test isperm(p) && isperm(q)
  @test v[p] == v[sortperm(v)] && v[q] == v[sortperm(v)]
end
