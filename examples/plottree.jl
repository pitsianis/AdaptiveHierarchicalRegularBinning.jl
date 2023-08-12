using PythonPlot, AdaptiveHierarchicalRegularBinning, AbstractTrees
import AdaptiveHierarchicalRegularBinning: center, sidelength, isleaf, children

function plottree(ax,node; dim1 = 1, dim2 = 2)
  c = center(node)
  h = sidelength(node) / 2
  x = [c[dim1]-h, c[dim1]+h, c[dim1]+h, c[dim1]-h, c[dim1]-h]
  y = [c[dim2]-h, c[dim2]-h, c[dim2]+h, c[dim2]+h, c[dim2]-h]

  if isleaf(node)
    ax.plot(x, y, color="red", linewidth=0.4)
  else
    # ax.plot(x, y, color="blue", linewidth=0.4)
    for child in children(node)
      plottree(ax, child; dim1, dim2)
    end
  end
  nothing
end

function plotpointstree(tree)
  fig, axs = subplots(layout="constrained", figsize=(10,10))
  axs.cla(); plottree(axs, tree)
  X = points(tree)
  axs.scatter(X[1,:], X[2,:], s=0.5, color="black")
  display( fig )
end

