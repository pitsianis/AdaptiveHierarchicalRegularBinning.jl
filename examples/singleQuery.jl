## A snipet to show how to perform a single query
using AdaptiveHierarchicalRegularBinning, AbstractTrees, Distances
±(x, y) = (x+y, x-y)
±(dx) = (x) -> x ± dx
square(x) = x .^ 2

function point2boxdist(p, node)
  h = sidelength(node) / 2   # Get the half-side lengths of the box
  c = center(node) # Get the center coordinates of the box

  # Calculate the vector from the point to the box's center
  v = abs.(p .- c)

  # Calculate the minimum distance for each component of the vector
  min_distance_vector = max.(v .- h, 0.0)

  # Calculate the magnitude of the minimum distance vector
  min_distance = sqrt(sum(min_distance_vector .^ 2))

  return min_distance
end

function predicate(node, query, r)
  return r >= point2boxdist(query, node)
end


function searchtree(tree, query, r=Inf, i=-1)
  if isleaf(tree)
    (d,k) = findmin( x -> evaluate(Euclidean(), x, query), eachcol( points(tree) ) )
    d < r && return ( d, range(tree)[k] )
  else
    for child in children(tree; by = x -> point2boxdist(query, x) )
      predicate(child, query, r) || continue
      r, i = searchtree(child, query, r, i)
    end
  end
  return (r, i)
end

using Test
function test()

  @testset "Boundary test" begin
      # TODO
  end

  @testset "Rand test" begin
    @testset "Iteration $i" for i in 1:10
      d = 6; n = 40_000
      X = randn(d, n)
      q = randn(d, 1) / 2rand()

      println("AHRB")
      @time begin
        tree = @time ahrb(X, 8, 1000; dims=2, QT=UInt)
        r, j = @time searchtree(tree, q, Inf, -1)
        i = tree.info.perm[j]
      end

      println("Brute force")
      result = j
      expected = @time argmin(reshape(sqrt.(sum((points(tree) .- q).^2,dims=leaddim(tree) == 1 ? 2 : 1)), :))
      println("result=$(result|>Int), expected=$(expected|>Int)")
      @test result == expected
    end
  end
end

test()
