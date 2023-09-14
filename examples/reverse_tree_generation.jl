using AbstractTrees
using StatsBase
using Random

mutable struct Node
  level
  code

  points

  children
end

Node(l) = Node(l, -1, nothing, Node[])


"""
    build_tree(l, d, bpl)

Builds a hierarchical tree of nodes.
This function recursively constructs the tree by distributing boxes (nodes) among levels, ensuring that each level can contain the specified number of boxes.

Parameters:
- `l::Int`: The current level of the tree.
- `d::Int`: The number of dimensions.
- `bpl::Vector{Int}`: An array specifying the number of boxes at each level of the tree.

Returns:
- `Vector{Node}`: An array of nodes representing the hierarchical tree.
"""
function build_tree(l, d, bpl)
  nb = bpl[l]
  L = length(bpl)
  if l == L
    nb > BigInt(2)^(d*l) && throw("Not enough boxes on level $l to contain $(bpl[l]) boxes")

    return [Node(l) for _ in 1:nb]
  else
    level_children = build_tree(l+1, d, bpl)
    non_terminal_nodes = []

    while !isempty(level_children)
      can_add_to = enumerate(non_terminal_nodes) |>
                      (X) -> Iterators.filter(x->length(last(x).children) < BigInt(2)^d, X) |>
                      (X) -> map(first, X)

      if length(non_terminal_nodes) == BigInt(2)^(l*d)-nb
        if isempty(can_add_to)
          throw("Not enough boxes on level $l to contain $(bpl[l]) boxes")
        end
      else
        push!(can_add_to, length(non_terminal_nodes)+1)
      end

      to_add_to = rand(can_add_to)

      if to_add_to == length(non_terminal_nodes)+1
        push!(non_terminal_nodes, Node(l))
      end

      free_space = BigInt(2)^d - length(non_terminal_nodes[to_add_to].children)
      to_add = min(rand(1:free_space), length(level_children))

      if isempty(non_terminal_nodes[to_add_to].children)
        to_add = min(rand(2:free_space), length(level_children))
        if to_add == 1
          push!(non_terminal_nodes[to_add_to].children, Node(l+1))
        end
      else
        to_add = min(rand(1:free_space), length(level_children))
      end
      append!(non_terminal_nodes[to_add_to].children, level_children[1:to_add])
      deleteat!(level_children, 1:to_add)
    end


    @assert all(2 <= length(n.children) <= 2^d for n in non_terminal_nodes)

    terminal_nodes = [Node(l) for _ in 1:nb]
    ret = [terminal_nodes..., non_terminal_nodes...]
    length(ret) > BigInt(2)^(l*d) && throw("Not enough boxes on level $l to contain $(bpl[l]) boxes")
    return ret
  end
end


"""
    assign_codes(tree::Node, d)

Assigns Morton codes to nodes in the tree.
This function assigns unique Morton codes to nodes in the tree.

Parameters:
- `tree::Node`: The root node of the hierarchical tree.
- `d::Int`: The number of dimensions.
"""
function assign_codes(tree::Node, d)

  length(tree.children) == 0 && return

  @assert length(tree.children) >= 2
  codes = zeros(Int, length(tree.children))
  codes[end] = 2^d-1
  if length(tree.children) > 2
    codes[2:end-1] = sample(1:2^d-2, length(tree.children)-2;replace=false)
  end


  for (code, node) in zip(codes, tree.children)
    node.code = code
    assign_codes(node, d)
  end

end



"""
    populate_tree(tree::Node, d, np, L, acc_codes=tuple())

Populates the tree with points.
This function populates the nodes in the tree with random points according to their Morton codes.

Parameters:
- `tree::Node`: The root node of the hierarchical tree.
- `d::Int`: The number of dimensions.
- `np::Int`: The maximum number of points per node to generate. (Minimum is np/2)
- `L::Int`: The number of levels in the tree.
- `acc_codes::Tuple{Int}`: A tuple containing accumulated node codes (internal use).
"""
function populate_tree(tree::Node, d, np, L, acc_codes=tuple())
  if isempty(tree.children)
    return populate_leaf(tree, d, np, acc_codes)
  end


  for node in tree.children
    populate_tree(node, d, np, L, (acc_codes..., node.code))
  end
end


function populate_leaf(tree::Node, d, np, codes)
  sidelen = 1.0 / (1 << length(codes))
  center = ones(d)
  for code in Iterators.reverse(codes)
    center /= 2
    for k in 1:d
        center[k] += (code >> (d - k)) & 1
    end
  end
  center /= 2

  # Generate np random points within the cube
  tree.points = rand(d, ceil(Int, max( cld(np,2)+1, rand() * np))) .* sidelen .+ center .- (sidelen/2)

  if all(codes .== 0)
    tree.points[:, 1] .= 0
  end

  if all(codes .== 2^d-1)
    tree.points[:, 1] .= 1
  end
  return
end

numpoints(tree) = isempty(children(tree)) ? size(tree.points, 2) : sum(numpoints, children(tree))
test_pmax(tree, np) =all((isempty(n.children) ? (numpoints(n)<= np) : (numpoints(n)>np))  for n in PostOrderDFS(tree))
test_unique(tree) = begin

  codes(node) = [n.code for n in children(node)]

  for node in PostOrderDFS(tree)
    c = codes(node)
    u = unique(c)
    if length(c) != length(u)
      return false
    end
  end
  return true
end
test_2branch(tree) = all( length(n.children) == 0 || length(n.children) >= 2  for n in PostOrderDFS(tree) )
test_contains_extrema(X) = any(all(col .== ones(length(col))) for col in eachcol(X)) && any(all(col .== zeros(length(col))) for col in eachcol(X))

AbstractTrees.children(tree::Node) = tree.children


function test_equal_trees(tree1, tree2)

  tree1.code == tree2.code || return false


  length(children(tree1)) != length(children(tree2)) && return false

  for (child1, child2) in zip(children(tree1), children(tree2))
    test_equal_trees(child1, child2) || return false
  end

  return true

end


"""
    gen_data(d, bpl, np; seed=0)

Generates data for a hierarchical tree.
This function creates a hierarchical tree, assigns codes to nodes, and populates the nodes with random points. It also performs several tests to ensure data integrity.

Parameters:
- `d::Int`: The number of dimensions.
- `bpl::Vector{Int}`: An array specifying the number of boxes at each level of the tree.
- `np::Int`: The maximum number of points per node to generate.
- `seed::Int`: (Optional) The random seed for reproducibility.

Returns:
- `Node`: The root node of the hierarchical tree.
- `Array{Float64, 2}`: An array of points generated for the tree.
"""
function gen_data(d, bpl, np; seed=0)

  Random.seed!(seed)

  tree = Node(0, 0, nothing, build_tree(1, d, bpl))
  assign_codes(tree, d)
  populate_tree(tree, d, np, length(bpl))

  X = mapreduce(t->t.points, hcat, Leaves(tree))

  @assert test_pmax(tree, np)
  @assert test_unique(tree)
  @assert test_2branch(tree)
  @assert test_contains_extrema(X)
  return tree, X
end