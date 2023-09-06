
using AdaptiveHierarchicalRegularBinning, AbstractTrees, Random, Profile

mutable struct NNinfo
  adjacentlist::Vector{Int64}
end
NNinfo() = NNinfo(Int64[])

function test()
  Random.seed!(0)
  n = 100_000
  d = 4
  maxL = min(120 ÷ d, 25)
  maxP = Int(ceil(sqrt(n)))

  X = rand(d, n)
  tree = AdaptiveHierarchicalRegularBinning.ahrb(X, maxL, maxP; QT=UInt128, ctxtype=NNinfo)
  println(tree)

  # initialize node context
  tree.info.context .= NNinfo.()

  # set adjacentlist

  @time Base.Threads.@threads for t in collect(PreOrderDFS(tree))
    for s in PreOrderDFS(x -> qbox2boxdistInf(t, x) == 0, tree)
      if nindex(t) != nindex(s) &&             # not coincident and
         (qbox2boxdistInf(t, s) == 0 &&         # adjacent
          (depth(t) == depth(s) ||              # at same level or
           (depth(t) > depth(s) && isleaf(s)))) # s is a leaf and at a coarser level
        push!(getcontext(t).adjacentlist, nindex(s))
      end
    end
  end

end

test()
# Profile.clear_malloc_data()
# test()
