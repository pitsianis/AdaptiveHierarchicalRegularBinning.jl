
orderchildrenindices(node,refnode) = cindices(node)[sortperm(children(node), by = x -> qbox2boxdistInf(x, refnode))]

function multilevelinteractions(t,s, prunepredicate, processleafpair!, postconsolidate!; orderchildren = false)
  function dtt(t, s)

    if prunepredicate(t, s)
      return
    end
    if nindex(t) == nindex(s) # coincident
      if isleaf(t)                      # both are leaves
        processleafpair!(t, s)
        postconsolidate!(t)
      else
        @threading for c in cindices(t)            # interact with selves first
          dtt(SpatialTree(TreeInfo(t), c), SpatialTree(TreeInfo(t), c))
        end
        @threading for t_child in cindices(t)      # interact with others next
          for s_child in cindices(s)
            t_child != s_child && dtt(SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child))
          end
        end
      end
    else
      if isleaf(t) && isleaf(s)         # both are leaves, we are done
        processleafpair!(t, s)
        postconsolidate!(t)
      elseif isleaf(t)
        sci = orderchildren ? orderchildrenindices(s,t) : cindices(s)
        for s_child in sci
          dtt(t, SpatialTree(TreeInfo(s), s_child))
        end
      elseif isleaf(s)
        tci = orderchildren ? orderchildrenindices(t,s) : cindices(t)
        for t_child in tci
          dtt(SpatialTree(TreeInfo(t), t_child), s)
        end
      else
        tci, sci = orderchildren ? (orderchildrenindices(t,s), orderchildrenindices(s,t)) : (cindices(t), cindices(s))
        for t_child in tci
          for s_child in sci
            dtt(SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child))
          end
        end
      end
    end
  end

  dtt(t, s)
end

function prioritymultilevelinteractions(t, s, nodedist, prunepredicate, processleafpair!, postconsolidate!)
  
  pq = PriorityQueue{Tuple{SpatialTree,SpatialTree},Tuple{Float64,Int}}()
  counter = 0
  pq[t,s] = (nodedist(t,s), counter -= 1)
  while !isempty(pq)
    t,s = dequeue!(pq)
    if prunepredicate(t, s)
      # do nothing
    elseif nindex(t) == nindex(s) # coincident
      if isleaf(t)                      # both are leaves
        processleafpair!(t, s)
        postconsolidate!(t)
      else
        for t_child in cindices(t)      # interact with others next
          for s_child in cindices(s)
            # prioritymultilevelinteractions(SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child), nodedist, prunepredicate, processleafpair!, postconsolidate!)
            pq[SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child)] = (0.0, counter -= 1)
          end
        end
      end
    else
      if isleaf(t) && isleaf(s)         # both are leaves, we are done
        processleafpair!(t, s)
        postconsolidate!(t)
      elseif isleaf(t)
        for s_child in cindices(s)
          pq[t, SpatialTree(TreeInfo(s), s_child)] = (nodedist(t, SpatialTree(TreeInfo(s), s_child)), counter -= 1)
        end
      elseif isleaf(s)
        for t_child in cindices(t)
          pq[SpatialTree(TreeInfo(t), t_child), s] = (nodedist(SpatialTree(TreeInfo(t), t_child), s), counter -= 1)
        end
      else
        for t_child in cindices(t)
          for s_child in cindices(s)
            pq[SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child)] = 
              (nodedist(SpatialTree(TreeInfo(t), t_child), SpatialTree(TreeInfo(s), s_child)), counter -= 1)
          end
        end
      end
    end
  end
end
