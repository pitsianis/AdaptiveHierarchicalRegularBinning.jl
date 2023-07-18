using AdaptiveHierarchicalRegularBinning

X = rand(10_000, 2)
X = vcat(X,X[1:100:901,:])

@time u = hcat(unique(eachrow(X))...)'

# OK solve this using the tree

@time tree = regural_bin(UInt32, X, 8, 100; dims=1)