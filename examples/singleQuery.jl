## A snipet to show how to perform a single query
using AdaptiveHierarchicalRegularBinning, AbstractTrees, Distances

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
      d = 6; n = 400_000
      X = randn(d, n)
      q = randn(d, 1) / 2rand()

      println("AHRB")
      # run once to precompile
      tree = ahrb(X, 8, 1000; QT=UInt)
      @time begin
        tree = @time ahrb(X, 8, 1000; QT=UInt)
        r, j = @time searchtree(tree, q, Inf, -1)
        i = tree.info.perm[j]
      end

      println("Brute force")
      result = j
      expected = @time argmin(reshape(sqrt.(sum((points(tree) .- q).^2,dims=1)), :))
      println("result=$(result|>Int), expected=$(expected|>Int)")
      @test result == expected
    end
  end
end

test()
