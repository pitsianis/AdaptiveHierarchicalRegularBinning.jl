using AbstractTrees

struct Tree
  data_vec
  lo
  hi
  tree_vec
  idx
  children
  depth
end


Tree(root::Tree, idx) = Tree(root.data_vec,
                             root.tree_vec[idx].lo,
                             root.tree_vec[idx].hi,
                             root.tree_vec,
                             idx,
                             root.tree_vec[idx].children,
                             root.tree_vec[idx].depth)

AbstractTrees.children(t::Tree) = Tree.(Ref(t), t.children)
AbstractTrees.nodevalue(t::Tree) = t.data_vec[t.lo]


function make_tree(R)
  length_unique(R)
end

function length_unique(R)
  acc = zero(length(R))

  length(R) == 0 && return acc

  p = first(R)
  acc += one(acc)

  # TODO: Try binary search
  for r in R
    if r != p
      p = r
      acc += one(acc)
    end
  end

  return acc
end



function make_tree_impl(root::Tree, L)

  root.depth > L     && return
  root.lo == root.hi && return

  tag(x) = (x >> (8*sizeof(x) - root.depth*8)) & 0xFF

  lo = root.lo
  idx = root.idx

  while lo <= root.hi
    lo_tag = tag(root.data_vec[lo])

    hi = lo+1
    while hi <= root.hi && tag(root.data_vec[hi]) == lo_tag
      hi += 1
    end
    hi -= 1

    idx = length(root.tree_vec)+1
    child = Tree(root.data_vec, lo, hi, root.tree_vec, idx, [], root.depth+1)
    push!(root.children, idx)
    push!(root.tree_vec, child)

    make_tree_impl(child, L)

    lo = hi+1
  end

end