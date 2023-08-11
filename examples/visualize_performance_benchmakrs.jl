using DrWatson
@quickactivate "AHRB"

using AdaptiveHierarchicalRegularBinning, BenchmarkTools, Random, DataFrames, PythonPlot,
  ColorSchemes, Colors

str_gitcommit = "00436a2ccfe3a17"
str_gitcommit = "18befd3c5"

df = collect_results!( datadir("experiments"); black_list = [] )

# sequential
df_filter = filter( 
  [:distribution,:enctype,:gitcommit,:np] => 
    (distr,enc,gitcommit,np) -> np == 1 && distr == "uniform" && enc == UInt128 && startswith( gitcommit, str_gitcommit ),
  df)

df_sub = df_filter[:, [:n,:d,:benchmark]]

df_sub[!,:benchmark] = map( x-> minimum(x.times) / 1e9, df_sub[!,:benchmark] )
# df_sub[!,:benchmark] = map( x-> x.allocs, df_sub[!,:benchmark] )

disallowmissing!(df_sub)

M = Matrix(df_sub)
M = sortslices(M;dims=1)

xx = M[:,1]
yy = M[:,2]
zz = M[:,3] ./ minimum(M[:,3])

# parallel
df_filter = filter( 
  [:distribution,:enctype,:gitcommit,:np] => 
    (distr,enc,gitcommit,np) -> np == 8 && distr == "uniform" && enc == UInt128 && startswith( gitcommit, str_gitcommit ),
  df)

df_sub = df_filter[:, [:n,:d,:benchmark]]

df_sub[!,:benchmark] = map( x-> minimum(x.times) / 1e9, df_sub[!,:benchmark] )
# df_sub[!,:benchmark] = map( x-> x.allocs, df_sub[!,:benchmark] )

disallowmissing!(df_sub)

Mpar = Matrix(df_sub)
Mpar = sortslices(Mpar;dims=1)

axesborder(border) = [border[1],border[2],border[3]-border[1],border[4]-border[2]]

pyplot.rcParams["text.usetex"] = true

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
z_label = ax.set_zlabel("normalized time", rotation=90)

# ax.set_title("Scaling with number of particles and dimensions")

ax.set_xticks(unique(xx) .+ 0.5*1e6, map( x-> "$(round(Int64,x/1e6))", unique(xx) ))
ax.set_yticks(unique(yy) .+ 0.5, map( x-> "$(round(Int64,x))", unique(yy) ))
ax.set_zticks(0:2:30)

ax.view_init(elev=20, azim=-150)

ax.plot(range(5,32,22), range(5,32,22) .* .8 .+ 4, zs=1.1e7, zdir="x", color="black", linestyle="--")
ax.plot(range(0.7e6,1.1e7,22), range(0.7,11.1,22) .* 2.55 .+ 1, zs=32 , zdir="y", color="black", linestyle="--")

ax.set_xlim(0.7e6, 1.1e7)
ax.set_ylim(5, 32)
ax.set_zlim(0, 30)
fig2

## parallel to sequential

@assert xx == Mpar[:,1]
@assert yy == Mpar[:,2]
zz = M[:,3] ./ Mpar[:,3]

fig3 = pyplot.figure( figsize=(6,5.5) )
border = [-0.1,0,1.05,1]
ax0 = fig3.add_axes(axesborder(border))
ax0.set_axis_off()

ax = fig3.add_axes(ax0.get_position(),projection="3d",zorder = 1)

colors = ( ColorSchemes.colorschemes[:RdPu_9] );
q = ( maximum(zz) - minimum(zz) ) / ( length(colors) - 1 );
clr = colors[ round.( Int, (zz .- minimum(zz)) ./ q ) .+ 1 ];

ax.bar3d(xx, yy, zeros(size(zz)), 0.5 * 1e6 * ones(size(xx)), 0.8 * ones(size(yy)), zz, 
  zsort="average", edgecolor = "black", linewidth=0.5, color = "#" .* hex.(clr), alpha=1.0)

ax.set_xlabel("number of particles \$\\times 10^6\$")
ax.set_ylabel("number of dimensions")
ax.zaxis.set_rotate_label(false)
z_label = ax.set_zlabel("parallel speedup (8 cores)", rotation=90)

# ax.set_title("Scaling with number of particles and dimensions")

ax.set_xticks(unique(xx) .+ 0.5*1e6, map( x-> "$(round(Int64,x/1e6))", unique(xx) ))
ax.set_yticks(unique(yy) .+ 0.5, map( x-> "$(round(Int64,x))", unique(yy) ))
ax.set_zticks(0:1:7)

ax.view_init(elev=20, azim=-150)

ax.set_xlim(0.7e6, 1.1e7)
ax.set_ylim(5, 32)
ax.set_zlim(0, 7)
fig3

## parallel only

@assert xx == Mpar[:,1]
@assert yy == Mpar[:,2]
zz = Mpar[:,3] ./ minimum(Mpar[:,3])

fig4 = pyplot.figure( figsize=(6,5.5) )
border = [-0.1,0,1.05,1]
ax0 = fig4.add_axes(axesborder(border))
ax0.set_axis_off()

ax = fig4.add_axes(ax0.get_position(),projection="3d",zorder = 1)

colors = ( ColorSchemes.colorschemes[:Blues_9] );
q = ( maximum(zz) - minimum(zz) ) / ( length(colors) - 1 );
clr = colors[ round.( Int, (zz .- minimum(zz)) ./ q ) .+ 1 ];

ax.bar3d(xx, yy, zeros(size(zz)), 0.5 * 1e6 * ones(size(xx)), 0.8 * ones(size(yy)), zz, 
  zsort="average", edgecolor = "black", linewidth=0.5, color = "#" .* hex.(clr), alpha=1.0)

ax.set_xlabel("number of particles \$\\times 10^6\$")
ax.set_ylabel("number of dimensions")
ax.zaxis.set_rotate_label(false)
z_label = ax.set_zlabel("normalized time (8 cores)", rotation=90)

# ax.set_title("Scaling with number of particles and dimensions")

ax.set_xticks(unique(xx) .+ 0.5*1e6, map( x-> "$(round(Int64,x/1e6))", unique(xx) ))
ax.set_yticks(unique(yy) .+ 0.5, map( x-> "$(round(Int64,x))", unique(yy) ))
ax.set_zticks(0:1:16)

ax.view_init(elev=20, azim=-150)

ax.set_xlim(0.7e6, 1.1e7)
ax.set_ylim(5, 32)
ax.set_zlim(0, 16)
fig4

## save figure

# fig2.savefig("/tmp/ahrb-performance.pdf",bbox_inches="tight")
# run(`pdfcrop /tmp/ahrb-performance.pdf /tmp/ahrb-performance.pdf`)

# fig3.savefig("/tmp/ahrb-speedup.pdf",bbox_inches="tight")
# run(`pdfcrop /tmp/ahrb-speedup.pdf /tmp/ahrb-speedup.pdf`)

# fig4.savefig("/tmp/ahrb-performance-parallel.pdf",bbox_inches="tight")
# run(`pdfcrop /tmp/ahrb-performance-parallel.pdf /tmp/ahrb-performance-parallel.pdf`)

# xx = reshape( M[:,1], 10, 10 )
# yy = reshape( M[:,2], 10, 10 )
# zz = reshape( M[:,3], 10, 10 ) ./ 1e9

# @assert all( xx[1,:]' .== xx )
# @assert all( yy[:,1]  .== yy )
# ax.imshow(zz)

# CS = ax.contour(xx, yy, zz, levels = 100)
# ax.clabel(CS, inline=true, fontsize=10)
# fig

error("Stop here!")

## OLD 2D BAR

fig, ax = pyplot.subplots(layout="constrained")

ax.cla()
width = 0.05  # the width of the bars
multiplier = 0
offset = width * multiplier

x = range(0,length=10)

for i in axes(yy,1)

  global offset = width * multiplier
  rects = ax.bar(x .+ offset, zz[i,:], width, label="d = $(Int64(yy[i,1]))", zorder = 3)
  # ax.bar_label(rects, padding=3)
  global multiplier += 1

end

ax.set_ylabel("time (sec)")
ax.set_xlabel("number of particles")
ax.set_title("Performance benchmark of AHRB")
ax.set_xticks(x .+ offset/2, map( x-> "$(round(Int64,x/1e6)) M", xx[1,:] ))
ax.set_yticks(0:10:100)
ax.legend(loc="upper left", ncols=3)
ax.set_ylim(0, 100)
ax.grid(true, which = "major", axis="y", zorder = -1, linestyle='-', linewidth=0.5)
ax.grid(true, which = "minor", axis="y", zorder = -1, linestyle='-', linewidth=0.2)
ax.minorticks_on()

ax.plot( x .+ offset/2, [7.8 .* (1:1:1e1);], color="black", linestyle="--", zorder = 2)

fig

## MATLAB

using MATLAB

@mput xx yy zz

mat"""
figure(1); clf

bar3(zz)
"""