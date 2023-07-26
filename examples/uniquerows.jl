using AdaptiveHierarchicalRegularBinning, AbstractTrees

X = rand(10, 1_000_000)
X[1, :] .= 1:size(X, 2)

X = hcat(X,X[:, 1:10:end])

@time expected = hcat(unique(eachcol(X))...)

# OK solve this using the tree

@time tree = ahrb!(UInt128, X, 8, 128; dims=2)

function m_unique(tree::SpatialTree)
  leaves = Leaves(tree)
  nleaves = sum((_)->1, leaves)
  R = Vector{Matrix}(undef, nleaves)
  for (i, leaf) in enumerate(Leaves(tree))
    p = (points(leaf))
    c = eachcol(p)
    u = unique(c)
    v = hcat(u...)
    R[i] = v
  end

  hcat(R...)
end

result = @time m_unique(tree)

I = sortperm(result[1, :])
result = result[:, I]

I = sortperm(expected[1, :])
expected = expected[:, I]

println("Test:  $(result == expected)")
