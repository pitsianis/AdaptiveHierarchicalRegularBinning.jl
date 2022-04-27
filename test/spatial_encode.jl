using AdaptiveHierarchicalRegularBinning: spatial_encode
using Test

@testset "basics" begin
  x = [0.0 0; 0 2; 2 0; 2 2]

  c,d,s = spatial_encode(x,1)

  @test d ≈ [0.0 0.0]
  @test s ≈ 0.5
  @test all(c .== [0; 2; 1; 3])
end
