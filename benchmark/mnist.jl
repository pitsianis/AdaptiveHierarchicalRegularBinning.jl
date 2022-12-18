# %% Imports
using AdaptiveHierarchicalRegularBinning
using LinearAlgebra
using MLDatasets
using BenchmarkTools

# %% Load data
X = MNIST(split=:train).features
X = reshape(X, prod(size(X)[1:end-1]), :)
X = transpose(X)

# %% Lower dimensionality
U, S, V = svd(X)
Np = 100
Up, Sp, Vp = U[:, 1:Np], S[1:Np], V[1:Np, 1:Np]

Xp = Up * Diagonal(Sp) * (Vp)

# %% Resample
np = size(Xp, 1)
n  = np*50
I  = (rand(UInt, n) .% np) .+ 1
Xr = Xp[I, :]

# %% Reshape
Xre = collect(transpose(Xr))

# %% Spatial encode
V = Xre

nd = 8
Vv = @view V[1:nd, :]
R = Vector{UInt}(undef, n)
l = 8

@benchmark AdaptiveHierarchicalRegularBinning.spatial_encode!(R, Vv, l; dims=Val(2), center=false)

# %% Radix Sort
Va = copy(V)
Ra = copy(R)
Ia = collect(UInt, 1:n)

Vb = similar(V)
Rb = similar(R)
Ib = Vector{UInt}(undef, n)

P = zeros(Bool, n)

rsd = AdaptiveHierarchicalRegularBinning.RadixSortDetails(nd, 1, n; dims=2, sml_th=1, dpt_th=l+1)

@benchmark(
  AdaptiveHierarchicalRegularBinning.radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, allocator),
  setup=begin
    Va .= V
    Ra .= R
    Ia .= 1:length(R)

    P .= false
    allocator = AdaptiveHierarchicalRegularBinning.Allocator(eltype(Ia))
  end
)
