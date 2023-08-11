using AdaptiveHierarchicalRegularBinning: fast_spatial_encode!
using Test

@testset "basics" begin
  x = [0.0 0.0;
       0.0 2.0;
       2.0 0.0;
       2.0 2.0]

  x = permutedims(x)

  r = Array{UInt}(undef, size(x, 2))
  d,s = fast_spatial_encode!(r, x, 1)

  @test d ≈ [0.0, 0.0]
  @test s ≈ 2.0
  @test all(r .== [0, 2, 1, 3])
end
