using AdaptiveHierarchicalRegularBinning
using Test, Random

include("datageneration.jl")

#
# Hook into Pkg.test so that tests from a single file can be run.  For example,
# to run only the invariance tests, use:
#
#   Pkg.test("AdaptiveHierarchicalRegularBinning", test_args=["invariants"])
#
# To run with 4 threads:
#
#   Pkg.test("AdaptiveHierarchicalRegularBinning", julia_args=["-t 4"])
#
# The above two commands can be combined. Also the name of the package is optional. 
# It will run the currently active package if not specified.
#
enabled_tests = lowercase.(ARGS)
function addtests(fname)
    key = lowercase(splitext(fname)[1])
    if isempty(enabled_tests) || key in enabled_tests
      Random.seed!(0)
      @testset "$(titlecase(fname)) tests" begin
        include("$fname.jl")
      end
    end
end

const files = (
  "dualtreetraversal",
  "boxdistance",
  "test_dense_nodes",
  "invariants",
  "newtree",
  "spatial_encode",
  "newtree",
  "knn"
)

@testset "AdaptiveHierarchicalRegularBinning.jl" begin
  @info "Running with $(Threads.nthreads()) threads"
  @info "Running $(AdaptiveHierarchicalRegularBinning.multithreading() ? "with" : "without") multithreading"
  map(files) do f
    addtests(f)
  end
end
