using DrWatson
@quickactivate "AHRB"

using AdaptiveHierarchicalRegularBinning, BenchmarkTools, Random, DataFrames, PythonPlot,
  ColorSchemes, Colors

str_gitcommit   = "00436a2ccfe3a17"
str_gitcommit_2 = "18befd3c5"

str_gitcommit_newversion = "050ec9fba8818"

df = collect_results!( datadir("experiments"); black_list = [] )

# sequential old
df_filter = filter( 
  [:distribution,:enctype,:gitcommit,:np] => 
    (distr,enc,gitcommit,np) -> np == 1 && distr == "uniform" && enc == UInt128 && ( startswith( gitcommit, str_gitcommit ) || startswith( gitcommit, str_gitcommit_2 ) ),
  df)

df_sub = df_filter[:, [:n,:d,:benchmark]]

df_sub[!,:benchmark] = map( x-> minimum(x.times) / 1e9, df_sub[!,:benchmark] )
# df_sub[!,:benchmark] = map( x-> x.allocs, df_sub[!,:benchmark] )

disallowmissing!(df_sub)


# sequential new
df_filter = filter( 
  [:distribution,:enctype,:gitcommit,:np] => 
    (distr,enc,gitcommit,np) -> np == 1 && distr == "uniform" && enc == UInt128 && startswith( gitcommit, str_gitcommit_newversion ),
  df)

df_sub_new = df_filter[:, [:n,:d,:benchmark]]

df_sub_new[!,:benchmark] = map( x-> minimum(x.times) / 1e9, df_sub_new[!,:benchmark] )
# df_sub[!,:benchmark] = map( x-> x.allocs, df_sub[!,:benchmark] )

disallowmissing!(df_sub)

## process

df_comp = innerjoin( df_sub, df_sub_new; on = [:n => :n, :d => :d], renamecols = "_old" => "_new" )

M = Matrix(df_comp)
M = sortslices(M;dims=1)

xx = M[:,1]
yy = M[:,2]

## figure
zz = M[:,3] ./ M[:,4]

axesborder(border) = [border[1],border[2],border[3]-border[1],border[4]-border[2]]
pyplot.rcParams["text.usetex"] = true

fig = pyplot.figure( figsize=(6,5.5) )
border = [-0.1,0,1.05,1]
ax0 = fig.add_axes(axesborder(border))
ax0.set_axis_off()

ax = fig.add_axes(ax0.get_position(),projection="3d",zorder = 1)

colors = ( ColorSchemes.colorschemes[:RdPu_9] );
q = ( maximum(zz) - minimum(zz) ) / ( length(colors) - 1 );
clr = colors[ round.( Int, (zz .- minimum(zz)) ./ q ) .+ 1 ];

ax.bar3d(xx, yy, zeros(size(zz)), 0.5 * 1e6 * ones(size(xx)), 0.8 * ones(size(yy)), zz, 
  zsort="average", edgecolor = "black", linewidth=0.5, color = "#" .* hex.(clr), alpha=1.0)

ax.set_xlabel("number of particles \$\\times 10^6\$")
ax.set_ylabel("number of dimensions")
ax.zaxis.set_rotate_label(false)
z_label = ax.set_zlabel("speedup \$\\frac{t_\\mathrm{old}}{t_\\mathrm{new}}\$", rotation=90)

# ax.set_title("Scaling with number of particles and dimensions")

ax.set_xticks(unique(xx) .+ 0.5*1e6, map( x-> "$(round(Int64,x/1e6))", unique(xx) ))
ax.set_yticks(unique(yy) .+ 1.5, map( x-> "$(round(Int64,x))", unique(yy) ))
ax.set_zticks(0:1:14)

ax.view_init(elev=20, azim=-150)

ax.set_xlim(0.7e6, 1.1e7)
ax.set_ylim(2, 32)
ax.set_zlim(0, 14)
fig

## figure just sequential new
zz = M[:,4]

fig2 = pyplot.figure( figsize=(6,5.5) )
border = [-0.1,0,1.05,1]
ax0 = fig2.add_axes(axesborder(border))
ax0.set_axis_off()

ax = fig2.add_axes(ax0.get_position(),projection="3d",zorder = 1)

colors = ( ColorSchemes.colorschemes[:Blues_9] );
q = ( maximum(zz) - minimum(zz) ) / ( length(colors) - 1 );
clr = colors[ round.( Int, (zz .- minimum(zz)) ./ q ) .+ 1 ];

ax.bar3d(xx, yy, zeros(size(zz)), 0.5 * 1e6 * ones(size(xx)), 0.8 * ones(size(yy)), zz, 
  zsort="average", edgecolor = "black", linewidth=0.5, color = "#" .* hex.(clr), alpha=1.0)

ax.set_xlabel("number of particles \$\\times 10^6\$")
ax.set_ylabel("number of dimensions")
ax.zaxis.set_rotate_label(false)
z_label = ax.set_zlabel("absolute time (seconds)", rotation=90)

# ax.set_title("Scaling with number of particles and dimensions")

ax.set_xticks(unique(xx) .+ 0.5*1e6, map( x-> "$(round(Int64,x/1e6))", unique(xx) ))
ax.set_yticks(unique(yy) .+ 1.0, map( x-> "$(round(Int64,x))", unique(yy) ))
ax.set_zticks(0:2:20)

ax.view_init(elev=20, azim=-150)

ax.plot(range(2,32,22), range(2,32,22) .* .25 .+ 12, zs=1.1e7, zdir="x", color="black", linestyle="--")
ax.plot(range(0.7e6,1.1e7,22), range(0.7,11.1,22) .* 1.80 .+ 0, zs=32 , zdir="y", color="black", linestyle="--")

ax.set_xlim(0.7e6, 1.1e7)
ax.set_ylim(2, 32)
ax.set_zlim(0, 20)
fig2

## export figures

# fig.savefig("/tmp/performance-comparison.png", bbox_inches="tight", dpi=300)
# fig2.savefig("/tmp/new-code-timing.png", bbox_inches="tight", dpi=300)
