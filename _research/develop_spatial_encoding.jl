using AdaptiveHierarchicalRegularBinning, BenchmarkTools, PythonPlot,
  ColorSchemes, Colors, Printf

import AdaptiveHierarchicalRegularBinning: spatial_encode!, fast_spatial_encode!

function benchmark_functions(n, d, maxL)
  
  V = rand(d, n); dims = 2

  ## test old function
  R = Vector{UInt128}(undef, size(V, dims))
  bitlen = size(V, dims==1 ? 2 : 1)
  maxL * bitlen > sizeof(UInt128) * 8 && throw("Not enough bits to represent the requested tree")

  ## test new function
  Rt = Vector{UInt128}(undef, size(V, dims))

  tslow = @benchmark spatial_encode!($R, $V, $maxL; dims=Val($dims), center=false);
  tfast = @benchmark fast_spatial_encode!($Rt, Vs, $maxL) setup=( Vs=copy($V) );

  @info maximum( abs.(R - Rt) )

  return [tslow, tfast]

end

T = Vector{BenchmarkTools.Trial}[]

for d = [2, 4, 8, 16, 32]
  maxL = 32÷d
  push!(T, benchmark_functions( 1_000_000, d, maxL ))
end

## plot

pyplot.rcParams["text.usetex"] = true
fig = pyplot.figure( figsize=(6,4) )
ax  = fig.gca()

ax.set_xscale("log", base=2)
ax.set_yscale("log", base=10)

ax.plot( [2, 4, 8, 16, 32], [minimum(t[1].times)/1e9 for t in T], "o-", label="spatial_encode!" )
ax.plot( [2, 4, 8, 16, 32], [minimum(t[2].times)/1e9 for t in T], "x-", label="fast_spatial_encode!" )

ax.set_xticks( [2, 4, 8, 16, 32],
  map( x -> @sprintf( "(%d, %d)", x, 32÷x ), [2, 4, 8, 16, 32] ) )
ax.set_xlabel("(number of dimensions, maximum tree level)")

ax.grid(true, which = "major", axis="y", zorder = -1, linestyle='-', linewidth=0.5)
ax.grid(true, which = "minor", axis="y", zorder = -1, linestyle='-', linewidth=0.2)

ax.set_ylabel("time [sec]")

for (i,d) = enumerate( [2, 4, 8, 16, 32] )
  maxL = 32÷d
  ytop = minimum(T[i][1].times)/1e9
  ybot = minimum(T[i][2].times)/1e9
  ymean = 10^( (log10(ytop) + log10(ybot))/2 )
  ax.text(d*1.02, ymean, @sprintf("x%.1f", ytop/ybot), horizontalalignment="left", verticalalignment="center")
  ax.plot( [d, d], [ybot, ytop], "k--", linewidth=0.5 )
end

ax.legend(loc="upper left")

# fig.savefig("/tmp/performance.png", bbox_inches="tight", dpi=300)

fig