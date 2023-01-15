# JuliaCon 2023 Proposal

## Abstract
`AdaptiveHierarchicalRegularBinning.jl` computes a hierarchical space-partitioning tree for a given set of points of arbitrary dimensions, that divides the space and stores the reordered points offering efficient access. Space-partitioning data structures are vital for algorithms that exploit spatial distance to reduce computational complexity, see for example the Fast Multipole Method, and algorithms finding nearest neighbors and their applications.


## Description
### Model
Assuming a set of `n` points `V` in a `d`-dimensional space, we partition the space into regular hierarchical bins. The resulting data structure is a sparse tree `T` with, at most, `2^d` child nodes per node and a maximum depth of `L`.

#### Normalization
The given set of points `V` is mapped into a `d`-dimensional unit hypercube via an affine transformation that scales and translates `V` resulting in `Vn`.

#### Binning
We split each dimension of the unit hypercube in half and recursively continue splitting each resulting hypercube in the same manner. Each partition is called a bin. The subdivision stops at a maximum depth `L` or when a bin contains `k` or fewer points. The recursive splitting process is recorded as a hierarchical tree data structure.

#### 1D Encoding
We use Morton encoding to map the `d`-dimensional set `Vn` to a one-dimensional space-filling curve. Each point in `Vn` is assigned an index in the reduced space resulting in `R`. Elements of `R` are bit-fields, each one consisting of `L` groups of `d` bits. These groups describe the position of the point in the corresponding level of the tree `T`. Points described by morton indices with equal most signifficant digits belong in the same bin. Thus, sorting `R` results in `Rs` which defines the sparse tree.

### Implementation
Cache-locality and parallel programming are some of the techniques we used, to make our implementation performant.

#### Cache locality
Cache locality offers fast memory access that greatly improves the performance of our algorithm.

- The set of points `V`, and by extension `Vn`, is defined as a `Matrix` of size `(d, n)`. The leading dimension describes the number of dimensions in the hyperspace resulting in an access pattern that is more friendly to the cache since all points of `V` have their corresponding coordinates densely packed in memory.

- Sorting `R` and `Vn` offers a memory layout that is cache-friendly. Since we access `V` through the `T` tree, we only access points that are in the same bin. `Rs` describes a bin of `T` with a contiguous block of memory, thus preserving cache-locality.

- `T` is not a linked-tree. `T` is a tree stored densely in memory as a `Vector{Node}`. Each `Node` of `T` is aware of their children and parent using their indices in this dense `Vector{Node}`.

#### Parallel Partial Sorting
The Morton curve bit-field denotes the tree node of each point. The partial sorting using the Most Significant Digit (MSD) radix-sort, places the points to the corresponding bins. Points that fall within the same leaf node do not get sorted. The radix-sort runs in parallel: the partition of digits is done with a parallel count-sort, and then each digit subset is processed independently in parallel.

#### Adaptive Tree
Empty bins, that is, nodes that do not contain any points, are not stored or referenced explicitly.


### References
- Sun, X., & Pitsianis, N. P. (2001). A Matrix Version of the Fast Multipole Method. In SIAM Review (Vol. 43, Issue 2, pp. 289–300). Society for Industrial & Applied Mathematics (SIAM). https://doi.org/10.1137/s0036144500370835

- Curtin, R., March, W., Ram, P., Anderson, D., Gray, A. & Isbell, C.. (2013). Tree-Independent Dual-Tree Algorithms. <i>Proceedings of the 30th International Conference on Machine Learning</i>, in <i>Proceedings of Machine Learning Research</i> 28(3):1435-1443 Available from https://proceedings.mlr.press/v28/curtin13.html.

- Greengard, L. (1990). The numerical solution of the n‐body problem. Computers in physics, 4(2), 142-152.

- Erdelyi, B. (2013). The fast multipole method for N-body problems. In AIP Conference Proceedings. ADVANCED ACCELERATOR CONCEPTS: 15th Advanced Accelerator Concepts Workshop. AIP. https://doi.org/10.1063/1.4773727

- Peterson, L. E. (2009). K-nearest neighbor. Scholarpedia, 4(2), 1883.

- Cho, M., Brand, D., Bordawekar, R., Finkler, U., Kulandaisamy, V., & Puri, R. (2015). PARADIS: An efficient parallel algorithm for in-place radix sort. Proceedings of the VLDB Endowment, 8(12), 1518-1529.

- Zagha, M., & Blelloch, G. E. (1991, August). Radix sort for vector multiprocessors. In Proceedings of the 1991 ACM/IEEE conference on Supercomputing (pp. 712-721).

- Morton, G. M. (1966). A computer oriented geodetic data base and a new technique in file sequencing.
