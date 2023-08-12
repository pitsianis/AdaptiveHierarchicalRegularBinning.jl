using PythonPlot, Random

include("plottree.jl")

Random.seed!(0)
d = 6; n = 800
X = zeros(d, n)
b = n ÷ 2
for i = 1:3
  jj = (i-1)*b÷2 .+ (1:b)
  println("$i $(jj[1]):$(jj[end]))")
  X[2i-1:2i,jj] = randn(2,length(jj)) .+ 4*i
end
M = rand(d, n) .> 0.50
X .+= 4*rand(d, n) .* M
# for _ = 1:100
# X .+= 0.04*rand(d, n)
tree = ahrb(X, 6, 2^3; dims=2, QT=UInt128)

fig, axs = subplots(d, d, layout="constrained", figsize=(10,10))
for i = 1:d, j = 1:d
  ax = axs[i-1,j-1]
  ax.cla(); plottree(ax, tree; dim1 = i, dim2 = j)
  ax.scatter( X[i,:], X[j,:], color="black", s=0.1)
  ax.set_aspect("equal")
  for tick in ax.xaxis.get_major_ticks()
    tick.tick1line.set_visible(false)
    tick.tick2line.set_visible(false)
    tick.label1.set_visible(false)
    tick.label2.set_visible(false)
  end
  for tick in ax.yaxis.get_major_ticks()
    tick.tick1line.set_visible(false)
    tick.tick2line.set_visible(false)
    tick.label1.set_visible(false)
    tick.label2.set_visible(false)
  end
end
display( fig )

# fig.savefig("/tmp/projections-2d.png", bbox_inches="tight", dpi = 300)