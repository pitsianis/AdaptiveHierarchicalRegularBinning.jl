#!/usr/bin/env sh
#SBATCH --mem=50G  # memory per node
#SBATCH -p compsci
#SBATCH --ntasks-per-node=5
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH -o %x-%A-%a.out
#=

echo julia --compiled-modules=no $(scontrol show job=$SLURM_JOBID | awk -F= '/Command=/{print $2}')
julia --compiled-modules=no $(scontrol show job=$SLURM_JOBID | awk -F= '/Command=/{print $2}')
exit
# =#



using DrWatson
@quickactivate "AHRB"

# detect if using SLURM
const IN_SLURM = "SLURM_JOBID" in keys(ENV)

# load packages
using Distributed
IN_SLURM && using ClusterManagers

# Here we create our parallel julia processes
if IN_SLURM
  cpus_per_task = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
  ENV["JULIA_NUM_THREADS"] = cpus_per_task
  pids = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
  print("\n")
else
  pids = nothing # pids = addprocs()   # uncomment this if you want to run in parallel mode even without Slurm
end

@everywhere using AdaptiveHierarchicalRegularBinning, DrWatson, BenchmarkTools, Random
@everywhere using Distributed
IN_SLURM && @everywhere using ClusterManagers

@everywhere include( "benchmarkutils.jl" )

force = true
path = datadir("experiments")
savecmd = exp -> savename(exp; ignores = ["cpu"], allowedtypes = (Real, String, Symbol, DataType))
general_args = Dict(
  :n => [(1:10) .* 1_000_000;], 
  :d => [1:10;],
  :seed => 0,
  :enctype => [UInt128, UInt64],
  :l => Derived( [:enctype, :d], (enctype, d) -> (enctype == UInt128 ? 128 ÷ d : 64 ÷ d ) ),
  :p => 2^7,
  :np => Threads.nthreads(),
  :hostname => gethostname(),
  :cpu => Sys.cpu_info()[1].model,
  :distribution => ["normal", "uniform"]
)

list_experiments = dict_list(general_args)


println(workers())
flush(stdout)

# Here we ask everyone to say hi!
# Output will appear in julia_in_parallel.output
if length( workers() ) > 1
  g = @sync @distributed (vcat) for w in workers()
      worker_id = myid()
      worker_host = gethostname()
      nth = Threads.nthreads()
      "Hello! I'm worker number $worker_id, with $nth threads and I reside on machine $worker_host. Nice to meet you!"
  end

  for i in g
    println(i)
  end
  flush(stdout)
end

println("Starting the experiments")
flush(stdout)

pmap( list_experiments ) do experiment
  data, _ = @produce_or_load(run_experiment, experiment, path; filename = savecmd, 
    tag = true, storepatch=true, force = force, 
    wsave_kwargs = Dict(:compress => true) )
  flush(stdout); flush(stderr)
end

!isnothing(pids) && rmprocs(pids)
println("procs removed")
flush(stdout)