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
X = rand(n,d)
L = 20
B = 32
tree = regural_bin(X, L, B; dims=1);

# What is tree? what are its fields?
# How do we access the tree nodes?
# How do we access the data array?

# be informal but precise
# we need to decide on the data structures

# We should have a recursive tree with the same structure 
# at each node except perhaps for the leaf nodes
# that means the call to `regular_bin` should return 
# multiple outputs

```


