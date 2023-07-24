# collect invariants from the code
# that is, identities that interelate the datastructure
# for any input
# these will be used to test the correctness of the code

using AdaptiveHierarchicalRegularBinning, AbstractTrees, LinearAlgebra

maxP = 3
maxL = 8
d = 8
n = 10_000
X = randn(2,400)
Xcopy = copy(X)
tree = regular_bin(UInt128, X, maxL, maxP; dims=2)

## Invariants

# Original points are permuted
@assert Xcopy[:,tree.info.perm] == tree.info.points
@assert X === tree.info.points

# all leaves have up to p points except the ones at the maxL level
@assert all(isleaf(node) ? 
              size(points(node),2) <= maxP : 
              true 
            for node in PreOrderDFS(tree) if depth(node) < maxL )

# all leaves are leaves
@assert all(isleaf.(Leaves(tree)))

# relationship of quantized and actual box centers and sides
@assert all(qbox(node) ≈ tree.info.scale * box(node) for node in PreOrderDFS(tree))

# each node represents a contiquous group of points, groups are ordered in preorder DFS
@assert all(minimum(low.(children(node))) == low(node) 
            for node in PreOrderDFS(tree) if !isleaf(node))

# users can add application-specific information to the tree

# root is only the first node ## BROKEN
isroot(tree)