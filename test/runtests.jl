using AdaptiveHierarchicalRegularBinning
using Test

const files = (
    "spatial_encode",
)

@testset "AdaptiveHierarchicalRegularBinning.jl" begin
    @testset "$(titlecase(f)) tests" for f in files
        include("$f.jl")
    end
end
