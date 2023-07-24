# collect invariants from the code
# that is, identities that interelate the datastructure
# for any input
# these will be used to test the correctness of the code

using AdaptiveHierarchicalRegularBinning, AbstractTrees

p = 25
maxL = 8
X = randn(2,400)
tree = regular_bin(UInt128, X, maxL, p; dims=2)

# Original points are permuted
X .-  tree.info.points[:,tree.info.perm] 

# all leaves have up to p points except the ones at the maxL level
@assert all( isleaf(node) ? length(points(node)) <= p : true for node in PreOrderDFS(tree) if depth(node) < maxL )

# I am still confused about the tree and the node, 
# how do I adress a node and how do I adress the subtree
# starting from this node?