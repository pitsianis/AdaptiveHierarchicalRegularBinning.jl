# %% Imports
using AdaptiveHierarchicalRegularBinning
using BenchmarkTools

# %% Control
d = 8
l = 8
n = 10_000_000
rsd = AdaptiveHierarchicalRegularBinning.RadixSortDetails(d, 1, n; dims=2, sml_th=1, dpt_th=l+1)


# %% Benchmark
@benchmark(
  AdaptiveHierarchicalRegularBinning.radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, allocator),
  setup=begin
    Va = rand(d, n)
    Ra = rand(UInt, n)
    Ia = collect(UInt, 1:n)

    Vb = similar(Va)
    Rb = similar(Ra)
    Ib = similar(Ia)

    P = falses(n)
    allocator = AdaptiveHierarchicalRegularBinning.Allocator(eltype(Ia))
  end
)