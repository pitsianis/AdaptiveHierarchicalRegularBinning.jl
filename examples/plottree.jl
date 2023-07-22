using PythonPlot, AdaptiveHierarchicalRegularBinning, AbstractTrees
import AdaptiveHierarchicalRegularBinning: center, box, isleaf, children

function plottree(ax,node; dim1 = 1, dim2 = 2)
  c = AdaptiveHierarchicalRegularBinning.center(node)
  h = box(node)# *.98
  x = [c[dim1]-h, c[dim1]+h, c[dim1]+h, c[dim1]-h, c[dim1]-h]
  y = [c[dim2]-h, c[dim2]-h, c[dim2]+h, c[dim2]+h, c[dim2]-h]

  if isleaf(node)
    ax.plot(x, y, color="red", linewidth=0.4)
  else
    # ax.plot(x, y, color="blue", linewidth=0.4)
    for child in AdaptiveHierarchicalRegularBinning.children(node)
      plottree(ax, child; dim1, dim2)
    end
  end
  fig
end


d = 4; n = 800
X = zeros(d, n)
for i = 1:d-1
  jj = 
  X[i:i+1,:] = randn(2,n) .+ 4*(d-1)
end
X = randn(d, n)
# for _ = 1:100
# X .+= 0.04*rand(d, n)
tree = regular_bin(UInt128, X, 6, 2^3; dims=2)

fig, axs = subplots(d, d, layout="constrained", figsize=(10,10))
for i = 1:d, j = 1:i-1
  ax = axs[i-1,j-1]
  ax.cla(); plottree(ax, tree; dim1 = i, dim2 = j)
  ax.scatter( X[i,:], X[j,:], color="black", s=0.1)
  ax.set_aspect("equal")
  ax.set_axis_off()
end
display( fig )
# sleep(0.1)
# end