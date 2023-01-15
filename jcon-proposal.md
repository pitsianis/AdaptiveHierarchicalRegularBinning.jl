# JuliaCon 2023 Proposal

## Abstract
Space-partitioning data structures are a vital part of the HPC ecosystem. They trade off space complexity and a small runtime overhead for an overall runtime complexity reduction, which can greatly reduce the amount of time needed for a computation. `AdaptiveHierarchicalRegularBinning.jl` offers a hierarchical space-partitioning tree that divides a space of arbitrary dimensions.


## Description
### Model
Assuming a set of `n` points `V` that defines a `d`-dimensional space, we offer the means to segregate that space into regular hierarchical bins. The resulting data structure is a hierarchical tree `T` with, at most, `2^d` child nodes per node and a maximum depth of `L`.

#### Normalization
We confine the given set of points `V` into a `d`-dimensional unit hypercube `Vn` via an affine transformation that scales and translates `V`.

#### Binning
We split the unit hypercube `Vn` in half for each dimension and recursively continue splitting each resulting hypercube in the same manner. Each hypercube is called a bin. The recursion process stops at a maximum depth of `L` or when a hypercube contains `k` or fewer points. This recursive splitting process eventually creates the hierarchical tree structure.

##### 1D Encoding
We use Morton encoding to map the `d`-dimensional set `Vn` to a one-dimensional space-filling curve `R`. Since a Morton curve is well-defined in a unit hypercube, we use the index of each point in the curve instead of the actual position in the reduced space. Each element in `R` is a bit-field divided into `L` groups of `d` bits, and each group describes the position of the corresponding point in the corresponding level of the tree `T`. Sorting `R` results in `Rs`, which is the initial state of the `T` tree. All elements of `Rs` that reside in the same node of `T` define a coherent block of the `Rs` array. These blocks are mutually exclusive from one another if the corresponding nodes are at the same level of the `T` tree. We, also, define `Vns` which is the permuted version of `Vn` that complies with `Rs`.

### Implementation
We took great care in implementing a performant yet generalized code that applies many of `Julia`'s best practices. Cache-locality, parallel programming, specialized functions, and minimizing allocations are some of the tools and techniques we used, to improve the runtime of the implementation.

#### Specialization and Generalization
Performance critical functions specialize over several parameters in order for the `Julia` compiler to produce performant code. The user is free to use a plethora of input types and leading dimensions for their data that comply with certain abstract types. This generalization comes at a probable runtime cost since not all data types and leading dimensions are as performant. We propose the most optimal input types for most cases.

#### Cache locality
Cache locality offers fast memory access that greatly improves the performance of our algorithm.

- The hyperspace `V`, and by extension both `Vn` and `Vns`, is defined as a `Matrix` of size `(d, n)`. The leading dimension describes the number of dimensions in the hyperspace resulting in an access pattern that is more friendly to the cache since all points of `V` have their corresponding coordinates densely packed in memory.

- Sorting `R` and `Vn` has several benefits, one of which is the memory access pattern that it offers. Since we only access `V` through the `T` tree, we only access points that are in the same node. `Rs`, and by extension `Vns`, describes a node of `T` with a contiguous block of memory, making operations on nodes cache-friendly.

- `T` is not a linked-tree. `T` is a tree stored densely in memory as a `Vecor{Node}`. Each `Node` of `T` is aware of their children and parent using their indices in this dense `Vector{Node}`.

#### Parallel Sorting
The nature of the Morton curve lends itself to parallelization schemes. As stated in previous sections each element in the `R` vector is an index in the Morton curve. This bit-field can be further analyzed since it contains information about the node location of the corresponding point. A Morton curve index can be divided into `L` groups of `d` bits, each group describes the node location on the appropriate level. This property can be exploited using an MSD radix-sort, which can be parallelized. In each recursion the sub-arrays to be sorted are densely packed in memory and mutually exclusive from one another, eliminating problems like cache locality, data races and data dependencies.

#### Sparse Tree
Assuming the depth of `T` is `L` and the tree is dense, then the tree is composed of `L` levels with `2^d` child nodes per node leading to a total of `(2^d)^L` nodes. This memory requirement is impractical for larger values of `d` and `L`. One way to reduce this memory requirement is for the `T` tree to be sparse. In a sparse tree, nodes that do not contain any points are not stored in memory, since we would never need to operate on such a node. A sparse tree has the same space complexity of `O((2^d)^L)` as a dense tree would, but in real-world applications with real-world data, the actual memory requirement is significantly less since most of the lower-level nodes would be empty and thus not stored.


### References
- Sun, X., & Pitsianis, N. P. (2001). A Matrix Version of the Fast Multipole Method. In SIAM Review (Vol. 43, Issue 2, pp. 289–300). Society for Industrial & Applied Mathematics (SIAM). https://doi.org/10.1137/s0036144500370835

- Curtin, R., March, W., Ram, P., Anderson, D., Gray, A. & Isbell, C.. (2013). Tree-Independent Dual-Tree Algorithms. <i>Proceedings of the 30th International Conference on Machine Learning</i>, in <i>Proceedings of Machine Learning Research</i> 28(3):1435-1443 Available from https://proceedings.mlr.press/v28/curtin13.html.

- Cho, M., Brand, D., Bordawekar, R., Finkler, U., Kulandaisamy, V., & Puri, R. (2015). PARADIS: An efficient parallel algorithm for in-place radix sort. Proceedings of the VLDB Endowment, 8(12), 1518-1529.

- Zagha, M., & Blelloch, G. E. (1991, August). Radix sort for vector multiprocessors. In Proceedings of the 1991 ACM/IEEE conference on Supercomputing (pp. 712-721).

- Morton, G. M. (1966). A computer oriented geodetic data base and a new technique in file sequencing.
