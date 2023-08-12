#!/usr/bin/env sh
#SBATCH --mem=50G  # memory per node
#SBATCH -p compsci
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -o %x-%A-%a.out
#=

echo julia --project=. -t$SLURM_CPUS_PER_TASK --compiled-modules=no $(scontrol show job=$SLURM_JOBID | awk -F= '/Command=/{print $2}')
julia --project=. -t$SLURM_CPUS_PER_TASK --compiled-modules=no $(scontrol show job=$SLURM_JOBID | awk -F= '/Command=/{print $2}')
exit
# =#

using DrWatson
@quickactivate "AHRB"

using AdaptiveHierarchicalRegularBinning, BenchmarkTools, Random, DataFrames, PythonPlot, Printf

# run new benchmarks or only load and show results?
run_new_benchmarks = false

# if running new benchmarks, force overwrite of existing results?
force = false

# where to store/load results
path = datadir("juliacon-predetermined-tree")

if run_new_benchmarks

  include( "benchmarkutils.jl" )
  include( "reverse_tree_generation.jl" )

  savecmd = exp -> savename(exp; ignores = ["cpu"], allowedtypes = (Real, String, Symbol, DataType))
  general_args = Dict(
    :p            => round.( Int, range( 2^7, 2^9, 7 ) ), # [2^7, 2^8, 2^9], # 2^10], 
    :n            => nothing,
    :d            => [5;],
    :seed         => 0,
    :enctype      => [UInt128],
    :l            => 8,
    :lstep        => 4,
    :method       => ["block-ecp", "fixed-length"],
    :np           => Threads.nthreads(),
    :hostname     => gethostname(),
    :cpu          => Sys.cpu_info()[1].model,
    :distribution => "juliacon-paper"
  )

  list_experiments = dict_list(general_args)

  # Here we ask everyone to say hi!
  # Output will appear in julia_in_parallel.output

  map( list_experiments ) do experiment
    data, _ = @produce_or_load(run_experiment, experiment, path; filename = savecmd, 
      tag = true, storepatch=true, force = force, 
      wsave_kwargs = Dict(:compress => true) )
  end

end

## collect results
df = collect_results!( path; black_list = [] )
disallowmissing!(df)
df[!,:n] = Int64.(df[!,:n])
df[!,:time] = map( x-> minimum(x.times) / 1e9, df[!,:benchmark] )

sort!( df, [:n, :p] )

@assert all( df.gitcommit[1] .== df.gitcommit )

## plot results
pyplot.rcParams["text.usetex"] = true

fig, ax = pyplot.subplots(layout="constrained")

colors = Dict( 
  "block-ecp" => "b",
  "fixed-length" => "r"
)

for (key, dfg) in pairs( groupby( df, :method ) )
  ax.plot( dfg.n, dfg.time, "-o", label = key[:method], color = colors[key[:method]] )
end

function xtickstr(x,y)

  str_x = if x < 1e6
    @sprintf( "%d \\!K", ceil(Int, x / 1e3) )
  else
    @sprintf( "%.1f \\!M", x / 1e6)
  end

  @sprintf "(%s, %d)" str_x y

end

ax.set_ylabel("time (sec)")
ax.set_xlabel("number of particles \\& leaf population \$(n, p_c)\$")
ax.grid(true, which = "major", axis="y", zorder = -1, linestyle='-', linewidth=0.5)
ax.grid(true, which = "minor", axis="y", zorder = -1, linestyle='-', linewidth=0.2)
ax.minorticks_on()
ax.tick_params(axis="x", which="minor", bottom=false)

nticks = sort( unique( [df.n df.p]; dims = 1 ); dims = 1 )
ax.set_xticks( nticks[:,1], map( xtickstr, nticks[:,1], nticks[:,2] ))

ax.legend(loc="upper left", shadow=true)

fig

# fig.savefig("/tmp/block-ecp-vs-fixed-length.pdf",bbox_inches="tight")