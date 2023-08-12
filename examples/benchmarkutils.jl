function run_experiment( experiment )

  if !haskey( experiment, :method )
    experiment[:method] = "block-ecp"
  end
  if !haskey( experiment, :lstep )
    experiment[:lstep] = 2
  end

  @unpack n,d,seed,distribution,enctype,p,l,method,lstep = experiment

  @assert enctype == UInt128

  result = copy( experiment )

  Random.seed!( seed )
  
  X = if distribution == "uniform"
    rand( d, n )
  elseif distribution == "normal"
    randn( d, n )
  elseif distribution == "juliacon-paper"
    _, X = gen_data( d, [ 2^3, 2^7, 2^9, 2^12, 2^7, 2^5, 2^2, 2^2 ], p ; seed = seed );
    n = size(X, 2); result[:n] = n
    X
  end

  b = @benchmark ahrb($X, $l, $p; QT=$enctype, method=$method, lstep=$lstep) evals=1 samples=5 seconds=5 setup=( GC.enable(false) ) teardown=( GC.enable(true); GC.gc() )

  result[:benchmark] = b

  return tostringdict( result )

end 