function run_experiment( experiment )

  @unpack n,d,seed,distribution,enctype,p,l = experiment

  result = copy( experiment )

  Random.seed!( seed )
  
  X = if distribution == "uniform"
    rand( d, n )
  else
    randn( d, n )
  end

  b = @benchmark regular_bin($enctype, $X, $l, $p; dims=2) evals=1 samples=5 seconds=5

  result[:benchmark] = b

  return tostringdict( result )

end 