# functions for computing distances between boxes and points

using AdaptiveHierarchicalRegularBinning, AbstractTrees
using Test

"""
    box2boxdist(c1, h1, c2, h2)

Compute the distance between two boxes.

c denotes the center of the box, h denotes the half-width of the box.

"""
function box2boxdist(c1, h1, c2, h2)
  # Path: examples/boxdistances.jl
  v = abs.(c1 .- c2)

  # Calculate the minimum distance for each component of the vector
  min_distance_vector = max.(v .- h1 .- h2, 0.0)

  # Calculate the magnitude of the minimum distance vector
  min_distance = sqrt(sum(min_distance_vector .^ 2))

  return min_distance
end

"""
    point2boxdist(p, c, h)

Compute the distance between a point and a box.
p denotes the point.
c denotes the center of the box, h denotes the half-width of the box.
"""
function point2boxdist(p, c, h)
  # Path: examples/boxdistances.jl
  v = abs.(p .- c)

  # Calculate the minimum distance for each component of the vector
  min_distance_vector = max.(v .- h, 0.0)

  # Calculate the magnitude of the minimum distance vector
  min_distance = sqrt(sum(min_distance_vector .^ 2))

  return min_distance
end

function foo(c1,h1,c2,h2)
  sqrt(sum(max.(c2 .- h2 .- (c1 .+ h1), c1 .- h1 .- (c2 .+ h2), 0.0).^2))
end

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

@testset "Two box-box distance test" begin

  n = 100_000; d = 3; L = 20; maxP = 2^11; X = randn(d,n); tree = ahrb!(X, L, maxP; QT = UInt)
  L = collect(Leaves(tree));
  nl = length(L)
  D1 = zeros(nl,nl)
  D2 = zeros(nl,nl)
  @time for i = 1:nl-1
    l1 = L[i]
    c1 = center(l1); h1 = sidelength(l1)/2;
    for j = i+1:nl
      l2 = L[j]
      c2 = center(l2); h2 = sidelength(l2)/2;
      D1[i,j] = foo(c1,h1,c2,h2); 
    end
  end

  @time for i = 1:nl-1
    l1 = L[i]
    c1 = center(l1); h1 = sidelength(l1)/2;
    for j = i+1:nl
      l2 = L[j]
      c2 = center(l2); h2 = sidelength(l2)/2;
      D2[i,j] = box2boxdist(c1, h1, c2, h2);
    end
  end

  println("$(maximum(abs.(D1 .- D2)))")
end