using AdaptiveHierarchicalRegularBinning
using Test

const files = (
  "spatial_encode",
  "countsort",
  "radixsort",
  "tree",
  "bit_interleave",
)

@testset "AdaptiveHierarchicalRegularBinning.jl" begin
    @testset "$(titlecase(f)) tests" for f in files
      println("Testing $f")
        include("$f.jl")
    end
end
