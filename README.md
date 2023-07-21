# AdaptiveHierarchicalRegularBinning

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pitsianis.github.io/AdaptiveHierarchicalRegularBinning.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pitsianis.github.io/AdaptiveHierarchicalRegularBinning.jl/dev)
[![Build Status](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml?query=branch%3Amain)

Given a cloud of points in a potentially high dimensional metric space, we partition the space into a hierarchy of regular bins, by subdividing the space adaptively up to a maximum L levels. A bin is not subdivided if it is at the maximum level `L` or if it contains up to `B` points.

*We need to decide on reasonable data structures for the input and output! Let's use this space to collect ideas and concerns and finalize the specification.*

Issues to consider:

 * How do we expect point coordinates? `n x d` or `d x n` or `dims=` or struct of arrays or array of structs?
 * What if the point cloud is an array of `struct` with extra information in addition to point coordinates, for example, point mass, point acceleration, etc
 * The output is a sparse tree, empty bins are not stored. See for example, adaptive quad or oct-tree in two and three dimensions respectively.
 * Each tree node provides pointers/indices to a data array where the *original* point records have been moved to be consecutive in memory.
 * The binning permutation of the points needs to be provided so that the points may be returned to their original positions.
 * All tree nodes (including the interior) need to point out to data structures that will hold auxiliary summary information (bin center coordinates and side lenght, level etc) or a pointer for specific applications (FMM polynomial expansion for example).
 * Our solution will handle arbitrary dimensional points but we may provide specialized implementations for 1-, 2- , 3- and 4d point clouds.

We provide (or use a ready-made) Morton spatial encoding function.

** Example **

```julia
using AdaptiveHierarchicalRegularBinning
n = 1000
d = 3
X = rand(d, n)
L = 20
B = 32
tree = regular_bin(UInt64, X, L, B; dims=2);

# What is tree? what are its fields?
# tree is just an index paired with a `TreeInfo` struct
# The index part is associated with the node pointed by the tree
#    and the relevant information (such as children_indices,
#    parent_index)
# The TreeInfo is the struct that holds all the global information
#    of the tree.
#    - Vector of nodes
#    - Vector of children
#
#    - Input points
#    - Vector of encoded points
#    - Resulting permutation after the reordering
#
#    - Spatial scale and offset
#    - Max depth and bin size
# How do we access the tree nodes?
# We can use the following which returns the tree nodes in their
# vector form, although this is considered esoteric operation
# and thus we do not offer an accessor function
tree.info.nodes
# Accessing the root node pointed by the tree has an accessor function
NodeInfo(tree)
# NodeInfo objects contain the following information:
# - Unit range to index the point matrix
# - Node depth
# - Parent index
# To traverse the tree we can utilize AbstractTrees.jl
using AbstractTrees
Leaves(tree)
PreOrderDFS(tree)
PostOrderDFS(tree)
# How do we access the data array?
# Again the data array is not considered as "exported" and so we do not offer
# accessors
tree.info.points # This is the reference to the input cloud of points
points(tree) # This is a view to the section of the points associated with the node

# To access data associated with a given node take a look at the following:
using AbstractTrees

for leaf in Leaves(tree)
  points(leaf) # This gets the points in the leaf node
  length(leaf) # This gets the point count of the leaf same as size(points(leaf), enum_dim)
  depth(leaf) # The depth of the leaf
  nodevalue(leaf) # The morton code of the leaf
  parent(leaf) # Gets the parent of the leaf
  children(leaf) # Gets an iterator to the children of the leaf
  center(leaf) # Computes the spatial center of the leaf
  box(leaf) # Computes the spatial boundary of the leaf
end

# how to see the nodevalue of the root node 
nodevalue(tree)
# What about the nodevalues of the children of the root node?

# the node values of leaves works as expected
nodevalues(Leaves(tree))
# but this one does not work
nodevalues(tree)
# why?

## Setup a context that is a struct
struct TestContext
  boxrange::Matrix{Float64}
  interactionlist::Vector{Int}
end

# initialize it
setcontext!(c1,TestContext(zeros(2,2),Int[]))
# modify it
push!(getcontext(c1).interactionlist, 2, 15)
```


