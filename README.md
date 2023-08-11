# AdaptiveHierarchicalRegularBinning

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Lifecycle:Maturing](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl)
<!-- [![status](??)](??) -->

The package
[AdaptiveHierarchicalRegularBinning.jl](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl),
or AHRB (pronounced as "arb", with a silent "h") for short, is a `Julia` package for binning data
point features.

The primary goal of AHRB  is to support and facilitate multi-resolution analysis of
particle-particle relations or interactions, especially for near-neighbor location or search at
multiple spatial scales or for far-neighbor filtering. The bins of AHRB are thereby chosen to be
$d$-dimensional cubes, extending the quad tree and oct tree to a $d$-dimensional hierarchical data
structure.


## Summary

In the base case, given a set of $n$ point particles, or feature vectors, in a $d$-dimensional
metric space, AHRB constructs a tree hierarchy that sorts the particles into nested bin nodes, up to
a cut-off level $\ell_{c}$.  The bin nodes at each level correspond to non-overlapping
$d$-dimensional cubes of the same size, each bin containing at least one particle, each non-leaf bin
containing more than $p_c$ particles.  The choice of the geometric shape and partition parameters
$\ell_{c}$ and $p_{c}$ serve the purpose of facilitating downstream tasks that involve accurate
multi-resolution analysis of particle-particle relationships.  AHRB offers additional
functionalities, especially for near-neighbor extraction or far-neighbor filtering at various
spatial scales.  When the feature dimension is low or modest, AHRB is competitive in time and space
complexities with other Julia packages for recursive binning of particles into nested cubes.
Distinctively, AHRB is capable of accommodating higher-dimensional data sets, without suffering from
high-order or exponential growth in memory usage with the increase in dimension. We demonstrate the
basic functionalities of AHRB and some extended ones, provide guaranteed time and space
complexities, and present sequential and parallel executions times on benchmarking datasets.


## Installation

```julia
] add https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl
```

## Examples

```julia
using AdaptiveHierarchicalRegularBinning, AbstractTrees
n = 100_000
d = 20
X = rand(d, n)

maxL = 6
maxP = 32
tree = ahrb(X, maxL, maxP; QT=UInt128);
```

### Properties & invariantes

```
# Original points are permuted
@assert X[:, tree.info.perm] == points(tree)

# all leaves have up to p points except the ones at the maxL level
@assert all(size(points(node), 2) <= maxP
  for node in PreOrderDFS(tree) if depth(node) < maxL && isleaf(node))

# all leaves are leaves
@assert all(isleaf.(Leaves(tree)))

# relationship of quantized and actual box centers and sides
@assert all(qbox(node) ≈ tree.info.scale * box(node) for node in PreOrderDFS(tree))

# each node represents a contiquous group of points, groups are ordered in preorder DFS
@assert all(minimum(low.(children(node))) == low(node) &&
          maximum(high.(children(node))) == high(node)
          for node in PreOrderDFS(tree) if !isleaf(node))
```
### Annotate each tree-node with mass & center of mass
```julia
masses = rand(n);

getmass(node::SpatialTree) = getcontext(node)[:mass]
getcom(node::SpatialTree)  = getcontext(node)[:com]

function populate_tree_ctx!(tree)
  foreach(PostOrderDFS(tree)) do node
    if isleaf(node)
      vm = masses[ tree.info.perm[range(node)] ]
      com  = points(node) * vm
      mass = sum( vm )    
    else
      com = zeros(size(pos,1)); mass = 0.0;
      for child in children(node)
        com .+= getcom(child) .* getmass(child)
        mass += getmass(child)
      end
    end
    com ./= mass
    setcontext!(node, (; com = com, mass = mass))
  end
end
```

## Applications

### Barnes-Hut N-body simulation

See [examples](examples/barneshut.jl)

https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/assets/4839092/face4a3d-0f17-49d6-9a2e-2377999bef4a

***
**AHRB at JuliaCon 2023**<br/>

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/Hkta_AEv5sA/0.jpg)](https://www.youtube.com/live/Hkta_AEv5sA?feature=share&t=15880)

  26th July 2023, 15:30–16:00 (US/Eastern)

***

<!--
## How to cite

If you use this software, please cite the following paper:

```bibtex
@inproceedings{floros2023ahrb,
    author = {Floros, Dimitris and Skourtis, Antonios and Pitsianis, Nikos and Sun, Xiaobai},
    doi = {xx},
    booktitle = {Proceedings of the JuliaCon Conferences},
    month = {??},
    title = {{Adaptive Hierarchical Regular Binning of Data Point Features}},
    year = {??}
}
```
-->
