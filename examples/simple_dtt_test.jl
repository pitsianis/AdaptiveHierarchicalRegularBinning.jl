using NearestNeighbors, AdaptiveHierarchicalRegularBinning, AbstractTrees, Random
using Test


function dtt(t, s)

  if nindex(t) == nindex(s) # coincident
    if isleaf(t)                      # both are leaves
      global D
      D[range(s),range(t)] .+= 1
    else
      Threads.@threads for c in cindices(t)            # interact with selves first
        dtt(SpatialTree(TreeInfo(t), c), SpatialTree(TreeInfo(t), c))
      end
      Threads.@threads for t_child in cindices(t)      # interact with others next
        for s_child in cindices(s)
          t_child != s_child && dtt(SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child))
        end
      end
    end
  else
    if isleaf(t) && isleaf(s)         # both are leaves, we are done
      global D
      D[range(s),range(t)] .+= 1
    elseif isleaf(t)
      for s_child in cindices(s)#; by = x -> qbox2boxdist(t, x))
        dtt(t, SpatialTree(TreeInfo(s), s_child))
      end
    elseif isleaf(s)
      for t_child in cindices(t)#; by = x -> qbox2boxdist(x, s))
        dtt(SpatialTree(TreeInfo(t), t_child), s)
      end
    else
      for t_child in cindices(t)
        for s_child in cindices(s)#; by = x -> qbox2boxdist(t_child, x) )
          dtt(SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child))
        end
      end
    end
  end

end

## 
V = rand(5, 1_000);
tree = ahrb(V,20,8;QT = UInt128, ctxtype=Vector{Tuple{Int64, Float64}});

D = zeros(Int, 1000, 1000)
dtt(tree, tree)

@assert D == ones(1000,1000)
