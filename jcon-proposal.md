# JuliaCon 2023 Proposal

## Abstract
Space-partitioning data structures are a vital part of the HPC ecosystem. They trade off space complexity and a small runtime overhead for an overall runtime complexity reduction, which can greatly reduce the amount of time needed for a computation. `AdaptiveHierarchicalRegularBinning.jl` offers a hierarchical space-partitioning tree that segregates a space of arbitrary dimensions.


## Description
### Model
Assuming a set of `n` points that defines a `d` dimentional space we offer the means to segregate that space into regular hierarchical bins. The resulting data structure is a hierarchical tree `T` with, at most, `2^d` nodes per level and maximum depth of `L`.

#### Normalization
We confine the given set of points `V` to a unit `d` dimentional hypercube by applying an affine transformation that scales and translates `V`. We call this normalized hypercube `Vn`.

#### Binning
We split the unit hypercube `Vn` in half for each dimension and recursively continue splitting each resulting hypercube in the same way. Each hypercube is called a bin. The recursion process stops at a maximum level `L` or when a hypercube has less than `k` points.

This splitting process eventualy creates a hierarchical tree structure with, at most, `2^d` nodes per tree level and a maximum depth of `L`.

##### 1D Encoding
We use morton encoding to map the `d` dimentional set `Vn` to an one dimentional space filling-curve `R`. Since a morton curve is well defined in a unit hypercube, we use the index of each point in the curve instead of an actual position in the reduced space. Each element in `R` is a bit field segregated into `L` groups of `d` bits, and each group describes the position of the corresponding point in the corresponding level of the tree `T`. Sorting `R` results in `Rs`, which is the initial state of `T`. In `Rs` all elements that reside in the same node of `T` define a coherent and mutialy exclusive, from other nodes in the same level, block of the `Rs` array. We, also, define `Vns` that is a permuted version of `Vn` that complies with `Rs`.

### Implementation
TODO

#### Cache locality
Cache locality offers fast memory access that greatly improves the performance of our algorithm.

- The hyperspace `V`, and by extension both `Vn` and `Vns`, is defined as a `Matrix` of size `(d, n)`, the leading dimension describes the number of dimensions in the hyperspace. This results in an access pattern that is more friendly to the cache, since all points in `V` have their corresponding coordinates densly packed in memory.

- Sorting `R` and `Vn` have several benefits, one of which is the memory access pattern that it offers. Since we only access `V` through the `T` tree, we only access points that are in the same node. `Rs`, and by extension `Vns`, describes a node of `T` with a contiguous block of memory, making operations on nodes cache friendly.

- `T` is not a linked-tree. `T` is a tree stored in densly packed `Vecor{Node}`. Each `Node` of `T` is aware of their children and parent using their indices in this dense `Vector{Node}`

#### Parallel Sorting
The nature of the morton curve lends itself to parallelization schemes. As stated in previous sections each element in the `R` matrix is an index in the morton curve. This bit field can be further analyzed since it contains information about the node location of the corresponding point. A morton curve index can be segregated into `L` groups of `d` bits, each group describing the node location on the corresponding level. This property can be exploited using an MSD radixsort, which can be parallelized since in each recursion the sub-arrays to be sorted are densly packed in memory and mutually exclusive from other sub-arrays, eliminating problems like data races and data dependencies.

