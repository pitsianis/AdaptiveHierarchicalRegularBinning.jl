# AdaptiveHierarchicalRegularBinning

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pitsianis.github.io/AdaptiveHierarchicalRegularBinning.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pitsianis.github.io/AdaptiveHierarchicalRegularBinning.jl/dev)
[![Build Status](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/actions/workflows/CI.yml?query=branch%3Amain)

Given a cloud of points in a potentially high dimensional metric space, we partition the points into a hierarchy of regular bins (hypercubes), by subdividing each hypercube adaptively up to a maximum L levels. A bin is not subdivided if it is at the maximum level `maxL` or contains up to `maxP` points.

*This package is under heavy development, and the API is not stable yet.*

** Example **

```julia
using AdaptiveHierarchicalRegularBinning
n = 100_000
d = 20
X = rand(d, n)
maxL = 6
maxP = 32
tree = ahrb!(X, maxL, maxP; dims=2, QT=UInt128);
```




