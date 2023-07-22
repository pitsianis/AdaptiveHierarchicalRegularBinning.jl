## A snipet to show how to perform a single query

query = struct
  coords = [0.5, 0.5]
  kthdist = Inf
end

function point2boxDist(p, node)
  h = box(node)
  c = center(node)
  if all(abs.(p .- c) .<= h)
    return 0
  else
    # compute the distance to the closest corner of the box (c +/- h)
  end

end

function predicate(node, query)
  return query.kthdist > point2boxDist(query.coords, node)
end


function searchTree(tree, query)
    if isleaf(tree)
      updatesearch!(query, points(tree), query))
    else
      cc = children(tree)
      dd = distances(cc, query)
      visitorder = sortperm(dd)
      for i in visitorder
        if predicate(cc[i], query)
          searchTree(cc[i], query)
        end
      end
    end
end