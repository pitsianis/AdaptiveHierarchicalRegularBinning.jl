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

# # detect if using SLURM
# const IN_SLURM = "SLURM_JOBID" in keys(ENV)

# # load packages
# using Distributed
# IN_SLURM && using ClusterManagers

# # Here we create our parallel julia processes
# if IN_SLURM
#   cpus_per_task = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
#   ENV["JULIA_NUM_THREADS"] = cpus_per_task
#   pids = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
#   print("\n")
# else
#   pids = nothing # pids = addprocs()   # uncomment this if you want to run in parallel mode even without Slurm
# end

# @everywhere using DrWatson

# @everywhere begin
#   @quickactivate "AHRB"
using AdaptiveHierarchicalRegularBinning, BenchmarkTools, Random
#   using Distributed
# end
# IN_SLURM && @everywhere using ClusterManagers

# @everywhere
include( "benchmarkutils.jl" )

force = true
path = datadir("experiments")
savecmd = exp -> savename("version-aug7-block-ecp-2023-dimitris", exp; ignores = ["cpu"], allowedtypes = (Real, String, Symbol, DataType))
general_args = Dict(
  :n => [(1:10) .* 1_000_000;], 
  :d => [2:4:42;],
  :seed => 0,
  :enctype => [UInt128],
  :l => 3,
  :lstep => 3,
  :p => 2^7,
  :np => Threads.nthreads(),
  :hostname => gethostname(),
  :cpu => Sys.cpu_info()[1].model,
  :distribution => ["uniform"]
)

list_experiments = dict_list(general_args)

# Here we ask everyone to say hi!
# Output will appear in julia_in_parallel.output

map( list_experiments ) do experiment
  data, _ = @produce_or_load(run_experiment, experiment, path; filename = savecmd, 
    tag = true, storepatch=true, force = force, 
    wsave_kwargs = Dict(:compress => true) )
end
