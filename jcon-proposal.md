# JuliaCon 2023 Proposal

## Abstract
Space-partitioning data structures are a vital part of the HPC ecosystem. They trade off space complexity and a small runtime overhead for an overall runtime complexity reduction, which can greatly reduce the amount of time needed for a computation. `AdaptiveHierarchicalRegularBinning.jl` offers a hierarchical space-partitioning tree that divides a space of arbitrary dimensions.


## Description
### Model
Assuming a set of `n` points that defines a `d` dimentional space we offer the means to segregate that space into regular hierarchical bins. The resulting data structure is a hierarchical tree `T` with, at most, `2^d` nodes per level and a maximum depth of `L`.

#### Normalization
We confine the given set of points `V` into a `d` dimentional unit hypercube via an affine transformation that scales and translates `V`. We call this normalized hypercube `Vn`.

#### Binning
We split the unit hypercube `Vn` in half for each dimension and recursively continue splitting each resulting hypercube in the same way. Each hypercube is called a bin. The recursion process stops at a maximum level `L` or when a hypercube contains `k` or less points.

This splitting process eventualy creates a hierarchical tree structure with, at most, `2^d` nodes per tree level and a maximum depth of `L`.

##### 1D Encoding
We use morton encoding to map the `d` dimentional set `Vn` to an one dimentional space-filling curve `R`. Since a morton curve is well defined in a unit hypercube, we use the index of each point in the curve instead of the actual position in the reduced space. Each element in `R` is a bit field divided into `L` groups of `d` bits, and each group describes the position of the corresponding point in the corresponding level of the tree `T`. Sorting `R` results in `Rs`, which is the initial state of the `T` tree. All elements of `Rs` that reside in the same node of `T` define a coherent block of the `Rs` array. These blocks are mutually exclusive from one another. We, also, define `Vns` that is a permuted version of `Vn` which complies with `Rs`.

### Implementation
We took great care in implementing a performant and yet generalized code that applies many of `Julia`'s best practices. Cache-locality, parallel programming, specialized functions and minimizing allocations are some of the tools and techniques we used, to improve the runtime of the implementation.

#### Specialization and Generalization
Performance critical functions specialize over several parameters in order for the `Julia` compiler to produce performant code. The user is free to use a plethora of input types and leading dimensions for their data that comply to certain abstract types. This generalization comes at a probable runtime cost since not all data types and leading dimensions are as performant. We propose the most optimal input types for most cases.

#### Cache locality
Cache locality offers fast memory access that greatly improves the performance of our algorithm.

- The hyperspace `V`, and by extension both `Vn` and `Vns`, is defined as a `Matrix` of size `(d, n)`, the leading dimension describes the number of dimensions in the hyperspace. This results in an access pattern that is more friendly to the cache, since all points of `V` have their corresponding coordinates densly packed in memory.

- Sorting `R` and `Vn` have several benefits, one of which is the memory access pattern that it offers. Since we only access `V` through the `T` tree, we only access points that are in the same node. `Rs`, and by extension `Vns`, describes a node of `T` with a contiguous block of memory, making operations on nodes cache friendly.

- `T` is not a linked-tree. `T` is a tree stored in a densly in memory as a `Vecor{Node}`. Each `Node` of `T` is aware of their children and parent using their indices in this dense `Vector{Node}`.

#### Parallel Sorting
The nature of the morton curve lends itself to parallelization schemes. As stated in previous sections each element in the `R` vector is an index in the morton curve. This bit field can be further analyzed since it contains information about the node location of the corresponding point. A morton curve index can be divided into `L` groups of `d` bits, each group describes the node location on the appropriate level. This property can be exploited using an MSD radixsort, which can be parallelized. In each recursion the sub-arrays to be sorted are densly packed in memory and mutually exclusive from one another, eliminating problems like cache locality, data races and data dependencies.
