using DrWatson
@quickactivate "AHRB"

using AdaptiveHierarchicalRegularBinning, BenchmarkTools, Random, DataFrames, PythonPlot,
  ColorSchemes, Colors

str_gitcommit_newversion = "050ec9fba8818"
str_gitcommit_newversion = "4386fa50c5d4a"

df = collect_results!( datadir("experiments"); black_list = [] )

# sequential new
df_filter = filter( 
  [:distribution,:enctype,:gitcommit,:np] => 
    (distr,enc,gitcommit,np) -> np == 1 && distr == "uniform" && enc == UInt128 && startswith( gitcommit, str_gitcommit_newversion ),
  df)

df_sub_new = df_filter[:, [:n,:d,:benchmark]]

df_sub_new[!,:benchmark] = map( x-> minimum(x.times) / 1e9, df_sub_new[!,:benchmark] )
# df_sub[!,:benchmark] = map( x-> x.allocs, df_sub[!,:benchmark] )

disallowmissing!(df_sub_new)

M = Matrix(df_sub_new)
M = sortslices(M;dims=1)

xx = M[:,1]
yy = M[:,2]
zz = M[:,3]

## figure
axesborder(border) = [border[1],border[2],border[3]-border[1],border[4]-border[2]]
pyplot.rcParams["text.usetex"] = true


fig = pyplot.figure( figsize=(6,5.5) )
border = [-0.1,0,1.05,1]
ax0 = fig.add_axes(axesborder(border))
ax0.set_axis_off()

ax = fig.add_axes(ax0.get_position(),projection="3d",zorder = 1)

colors = ( ColorSchemes.colorschemes[:Blues_9] );
q = ( maximum(zz) - minimum(zz) ) / ( length(colors) - 1 );
clr = colors[ round.( Int, (zz .- minimum(zz)) ./ q ) .+ 1 ];

ax.bar3d(xx, yy, zeros(size(zz)), 0.5 * 1e6 * ones(size(xx)), 2.0 * ones(size(yy)), zz, 
  zsort="average", edgecolor = "black", linewidth=0.5, color = "#" .* hex.(clr), alpha=1.0)

ax.set_xlabel("number of particles \$\\times 10^6\$")
ax.set_ylabel("number of dimensions")
ax.zaxis.set_rotate_label(false)
z_label = ax.set_zlabel("absolute time (seconds)", rotation=90)

# ax.set_title("Scaling with number of particles and dimensions")

ax.set_xticks(unique(xx) .+ 0.5*1e6, map( x-> "$(round(Int64,x/1e6))", unique(xx) ))
ax.set_yticks(unique(yy) .+ 1.5, map( x-> "$(round(Int64,x))", unique(yy) ))
ax.set_zticks(0:2:30)

ax.view_init(elev=20, azim=-150)

ax.plot(range(2,44,22), range(2,44,22) .* .58 .+ 5, zs=1.1e7, zdir="x", color="black", linestyle="--")
ax.plot(range(0.7e6,1.1e7,22), range(0.7,11.1,22) .* 2.70 .+ 0, zs=44 , zdir="y", color="black", linestyle="--")

ax.set_xlim(0.7e6, 1.1e7)
ax.set_ylim(2, 44)
ax.set_zlim(0, 30)
fig

## export figures

# fig.savefig("/tmp/performance-single-block.pdf", bbox_inches="tight")

