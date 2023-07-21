using DrWatson, BenchmarkTools, Random

@quickactivate "AHRB"

function run_experiment( experiment )

  @unpack n,d,seed = experiment

  result = copy( experiment )

  Random.seed!( seed )


end 

force = true
savecmd = exp -> savename(exp; ignores = [])
general_args = Dict(
  :n => [(1:10) .* 1e6;], 
  :d => [1:10;],
  :seed => 0,
  :hostname => gethostname(),
  :cpu => Sys.cpu_info()[1].model
)

list_experiments = dict_list(general_args)

path = datadir("experiments")

for experiment in list_experiments
  @show experiment
  @produce_or_load(run_experiment, experiment, path; filename = savecmd, 
    tag = true, storepatch=true, force = force, 
    wsave_kwargs = Dict(:compress => true) )
end