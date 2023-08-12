using AdaptiveHierarchicalRegularBinning, AbstractTrees
using Test

@testset "Elementary point-box distance test" begin

  @test point2boxdist([0; 0], [0.5;0.5], 0.5) == 0.0
  @test point2boxdist([0; 1], [0.5;0.5], 0.5) == 0.0
  @test point2boxdist([1; 0], [0.5;0.5], 0.5) == 0.0
  @test point2boxdist([1; 1], [0.5;0.5], 0.5) == 0.0

  @test point2boxdist([0; 0.5], [0.5;0.5], 0.5) == 0.0
  @test point2boxdist([0.5; 0], [0.5;0.5], 0.5) == 0.0

  @test point2boxdist([0; 0], [0.5;0.5], 0.25) == sqrt(2*0.25^2)
  @test point2boxdist([0; 1], [0.5;0.5], 0.25) == sqrt(2*0.25^2)
  @test point2boxdist([1; 0], [0.5;0.5], 0.25) == sqrt(2*0.25^2)
  @test point2boxdist([1; 1], [0.5;0.5], 0.25) == sqrt(2*0.25^2)

end

@testset "Elementary box-box distance test" begin

  @test box2boxdist([0; 0], 0.5, [1.0; 1.0], 0.5) == 0.0

  @test box2boxdist([0; 0], 0.25, [0.75;0.75], 0.25) == sqrt(2*0.25^2)
  @test box2boxdist([1; 0], 0.25, [0.25;0.75], 0.25) == sqrt(2*0.25^2)
  @test box2boxdist([0; 1], 0.25, [0.75;0.25], 0.25) == sqrt(2*0.25^2)

  @test box2boxdist([1; 0], 0.125, [0.25;0.75], 0.25) == sqrt(2*0.375^2)

end

