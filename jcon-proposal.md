# JuliaCon 2023 Proposal

## Abstract
`AdaptiveHierarchicalRegularBinning.jl` computes a hierarchical space-partitioning tree for a given set of points of arbitrary dimensions, that divides the space and stores the reordered points offering efficient access. Space-partitioning data structures are vital for algorithms that exploit spatial distance to reduce computational complexity, see for example the Fast Multipole Method, and algorithms finding nearest neighbors and their applications.  


## Description
### Model
Assuming a set of `n` points `V` in a `d`-dimensional space, we partition the space into regular hierarchical bins. The resulting data structure is a sparse tree `T` with, at most, `2^d` child nodes per node and a maximum depth of `L`. 

#### Normalization
The given set of points `V` is mapped into a `d`-dimensional unit hypercube `Vn` via an affine transformation that scales and translates `V`.

#### Binning
We split each dimension of the unit hypercube `Vn` in half and recursively continue splitting each resulting hypercube in the same manner. Each partition is called a bin. The subdivision stops at a maximum depth `L` or when a bin contains `k` or fewer points. The recursive splitting process is recorded as a hierarchical tree data structure.

##### 1D Encoding
The map of the `d`-dimensional point set `Vn` to a `d`-dimensional unit hypercube corresponds to a one-dimensional Morton space-filling curve `R`. Each element in `R` is a bit-field consisting of `L` groups of `d` bits. Each group describes the position of the point in the corresponding level of the tree `T`. Sorting of the point codes onto `R`, results to the order of the points in the space tree.


#### Cache locality
Cache locality offers fast memory access that greatly improves the performance of our algorithm.

- The hyperspace `V`, and by extension both `Vn` and `Vns`, is defined as a `Matrix` of size `(d, n)`. The leading dimension describes the number of dimensions in the hyperspace resulting in an access pattern that is more friendly to the cache since all points of `V` have their corresponding coordinates densely packed in memory.

- Sorting `R` and `Vn` has several benefits, one of which is the memory access pattern that it offers. Since we only access `V` through the `T` tree, we only access points that are in the same node. `Rs`, and by extension `Vns`, describes a node of `T` with a contiguous block of memory, making operations on nodes cache-friendly.

- `T` is not a linked-tree. `T` is a tree stored densely in memory as a `Vector{Node}`. Each `Node` of `T` is aware of their children and parent using their indices in this dense `Vector{Node}`.

#### Parallel Partial Sorting
The  Morton curve bit-field denotes the tree node of each point. The partial sorting using the Most Significant Digit (MSD) radix-sort, places the points to the corresponding bins. Points that fall within the same bin do not get sorted. The radix sort runs in parallel: the partition of digits is done with a parallel count sort, and then each digit subset is processed independently in parallel.

#### Adaptive Tree
Empty bins, that is, nodes that do not contain any points, are not stored or referenced explicitly.


### References
- Sun, X., & Pitsianis, N. P. (2001). A Matrix Version of the Fast Multipole Method. In SIAM Review (Vol. 43, Issue 2, pp. 289–300). Society for Industrial & Applied Mathematics (SIAM). https://doi.org/10.1137/s0036144500370835

- Curtin, R., March, W., Ram, P., Anderson, D., Gray, A. & Isbell, C.. (2013). Tree-Independent Dual-Tree Algorithms. <i>Proceedings of the 30th International Conference on Machine Learning</i>, in <i>Proceedings of Machine Learning Research</i> 28(3):1435-1443 Available from https://proceedings.mlr.press/v28/curtin13.html.

- Cho, M., Brand, D., Bordawekar, R., Finkler, U., Kulandaisamy, V., & Puri, R. (2015). PARADIS: An efficient parallel algorithm for in-place radix sort. Proceedings of the VLDB Endowment, 8(12), 1518-1529.

- Zagha, M., & Blelloch, G. E. (1991, August). Radix sort for vector multiprocessors. In Proceedings of the 1991 ACM/IEEE conference on Supercomputing (pp. 712-721).

- Morton, G. M. (1966). A computer oriented geodetic data base and a new technique in file sequencing.
