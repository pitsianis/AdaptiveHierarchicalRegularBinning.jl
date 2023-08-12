using AdaptiveHierarchicalRegularBinning, AbstractTrees

using BenchmarkTools
using JLD
using Random: MersenneTwister

import Printf: @sprintf
import AdaptiveHierarchicalRegularBinning: fast_spatial_encode!, countandpermute!, make_tree

include("generate_report.jl")

const EXTENSIVE_BENCHMARK = false

const SUITE = BenchmarkGroup()
SUITE["ahrb"] = BenchmarkGroup()
SUITE["encode"] = BenchmarkGroup()
SUITE["countandpermute"] = BenchmarkGroup()
SUITE["make tree"] = BenchmarkGroup()
SUITE["dual tree traversal"] = BenchmarkGroup()

for n_points in (EXTENSIVE_BENCHMARK ? (10^5, 10^6) : 10^6)
  for dim in (EXTENSIVE_BENCHMARK ? (1:6) : (1, 3, 5))
    for maxL in (EXTENSIVE_BENCHMARK ? (5, 10, 15, 20) : 20)
      for maxP in (EXTENSIVE_BENCHMARK ? (1, 16, 128) : 128)

        SUITE["ahrb"]["$dim × $n_points, maxL = $maxL, maxP = $maxP"] = 
          @benchmarkable ahrb(X, $maxL, $maxP; QT=UInt128) setup=(
            X = rand(MersenneTwister(1), $dim, $n_points);
            GC.enable(false)
          ) teardown=(
            GC.enable(true); 
            GC.gc()
          )
        
        SUITE["encode"]["$dim × $n_points, maxL = $maxL"] = 
          @benchmarkable fast_spatial_encode!(R, V, $maxL) setup=(
            V = rand(MersenneTwister(1), $dim, $n_points); 
            R = zeros(UInt128, $n_points);
            GC.enable(false)
          ) teardown=(
            GC.enable(true); 
            GC.gc()
          )

        SUITE["countandpermute"]["$dim × $n_points, maxL = $maxL"] =
          @benchmarkable countandpermute!(ix, Vp, Rp, V, R) setup=(
            V = rand(MersenneTwister(1), $dim, $n_points); 
            R = zeros(UInt128, $n_points);
            fast_spatial_encode!(R, V, $maxL);
            ix = [LinearIndices(R);];
            Vp = copy(V);
            Rp = copy(R);
            GC.enable(false)
          ) teardown=(
            GC.enable(true); 
            GC.gc()
          )

        SUITE["make tree"]["$dim × $n_points, maxL = $maxL, maxP = $maxP"] =
          @benchmarkable make_tree(Vp, Rp, ix, $maxL, $maxP, scale, offset; ctxtype = Bool) setup=(
            V = rand(MersenneTwister(1), $dim, $n_points); 
            R = zeros(UInt128, $n_points);
            (offset, scale) = fast_spatial_encode!(R, V, $maxL);
            ix = [LinearIndices(R);];
            Vp = copy(V);
            Rp = copy(R);
            countandpermute!(ix, Vp, Rp, V, R);
            GC.enable(false)
          ) teardown=(
            GC.enable(true); 
            GC.gc()
          )

      end
    end

    for maxL in (EXTENSIVE_BENCHMARK ? (5, 10, 15, 20) : 20)
      for maxP in [floor( Int, 20*sqrt( n_points ) );]

        SUITE["dual tree traversal"]["$dim × $n_points, maxL = $maxL, maxP = $maxP"] =
          @benchmarkable dualtreetraversal(T, T) setup=(
            V = rand(MersenneTwister(1), $dim, $n_points); 
            T = ahrb(V, $maxL, $maxP; QT=UInt128, ctxtype=Vector{Tuple{Int64, Float64}});
            foreach(node -> setcontext!(node, Tuple{Int64, Float64}[]), PreOrderDFS(T));
            nl = leafcount(T);
            foreach(node -> sizehint!(getcontext(node), nl), PreOrderDFS(T));
          )

      end
    end

  end
end

function run_benchmarks(name)
  paramspath = joinpath(dirname(@__FILE__), EXTENSIVE_BENCHMARK ? "params_extensive.jld" : "params.jld")
  if !isfile(paramspath)
      println("Tuning benchmarks...")
      tune!(SUITE)
      JLD.save(paramspath, "SUITE", params(SUITE))
  end
  loadparams!(SUITE, JLD.load(paramspath, "SUITE"), :evals, :samples)
  results = run(SUITE, verbose = true, seconds = 2)
  name = EXTENSIVE_BENCHMARK ? @sprintf("%s_extensive", name) : name
  JLD.save(joinpath(dirname(@__FILE__), name * ".jld"), "results", results)
end

function generate_report(v1, v2)
  v1_res = load(joinpath(dirname(@__FILE__), v1 * ".jld"), "results")
  v2_res = load(joinpath(dirname(@__FILE__), v2 * ".jld"), "results")
  open(joinpath(dirname(@__FILE__), "results_compare.md"), "w") do f
      printreport(f, judge(minimum(v1_res), minimum(v2_res)); iscomparisonjob = true)
  end
end

function generate_report(v1)
  v1_res = load(joinpath(dirname(@__FILE__), v1 * ".jld"), "results")
  open(joinpath(dirname(@__FILE__), "results_single.md"), "w") do f
      printreport(f, minimum(v1_res); iscomparisonjob = false)
  end
end

# run_benchmarks("baseline")
# generate_report("baseline") # generate a report with stats about a run
# run_benchmarks("new")
# generate_report("new", "baseline") # generate report comparing two runs
