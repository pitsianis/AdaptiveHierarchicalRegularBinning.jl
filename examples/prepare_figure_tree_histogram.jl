using DrWatson
@quickactivate "AHRB"

using AdaptiveHierarchicalRegularBinning, BenchmarkTools, Random, DataFrames, PythonPlot, Printf, StatsBase

include( "reverse_tree_generation.jl" )

d    = 5
p    = 5
seed = 42
_, X = gen_data( d, [ 2^3, 2^7, 2^9, 2^12, 2^7, 2^5, 2^2, 2^2 ], p ; seed = seed );
tree = ahrb(X, 8, p; QT=UInt128, lstep=4)


depths    = map( x->x.dpt, tree.info.nodes )
idx_leaf  = findall( isempty, tree.info.children )
idx_nonleaf = findall( !isempty, tree.info.children )
n_leaf    = fit(Histogram, depths[idx_leaf], 1:1.0:9.0).weights
n_nonleaf = fit(Histogram, depths[idx_nonleaf], 1:1.0:9.0).weights

n_leaf    = n_leaf ./ sum(n_leaf)
n_nonleaf = n_nonleaf ./ sum(n_nonleaf)

## plot results
pyplot.rcParams["text.usetex"] = true

fig, ax = pyplot.subplots(layout="constrained")

colors = Dict( 
  "leaf" => "#1f77b4",
  "non-leaf" => "#ff7f0e"
)

width = 0.25  # the width of the bars
multiplier = 0

offset = width * multiplier
x = 1:length(n_leaf)
ax.bar(x .+ offset, n_leaf, width, label="leaf", color = colors["leaf"], zorder = 3)

multiplier += 1
offset = width * multiplier
ax.bar(x .+ offset, n_nonleaf, width, label="non-leaf", color = colors["non-leaf"], zorder = 3)

ax.set_ylabel("percentage")
ax.set_xlabel("level")
ax.grid(true, which = "major", axis="y", zorder = -1, linestyle="-" , linewidth=0.5)
ax.grid(true, which = "minor", axis="y", zorder = -1, linestyle="--", linewidth=0.5, alpha = 0.5)
ax.minorticks_on()
ax.tick_params(axis="x", which="minor", bottom=false)
ax.set_yscale("log")

nticks = x
ax.set_xticks( nticks )

ax.legend(loc="upper left", shadow=true)

fig

# fig.savefig("/tmp/tree_hist_boxes.pdf",bbox_inches="tight")