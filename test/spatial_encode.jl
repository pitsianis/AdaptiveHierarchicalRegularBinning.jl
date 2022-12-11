using AdaptiveHierarchicalRegularBinning: spatial_encode!
using Test

@testset "basics" begin
  x = [0.0 0.0;
       0.0 2.0;
       2.0 0.0;
       2.0 2.0]

  r = Array{UInt}(undef, size(x, 1))
  d,s = spatial_encode!(r, x, 1; dims=Val(1), center=false)

  @test d ≈ [0.0, 0.0]
  @test s ≈ 0.5
  @test all(r .== [0, 2, 1, 3])
end
