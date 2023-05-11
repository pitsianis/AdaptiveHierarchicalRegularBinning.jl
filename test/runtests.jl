using AdaptiveHierarchicalRegularBinning
using Test

const files = (
  "bit_interleave",
  "spatial_encode",
  "countsort",
  "radixsort",
  "tree",
)

@testset "AdaptiveHierarchicalRegularBinning.jl" begin
    @testset "$(titlecase(f)) tests" for f in files
        include("$f.jl")
    end
end
