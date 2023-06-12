var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = AdaptiveHierarchicalRegularBinning","category":"page"},{"location":"#AdaptiveHierarchicalRegularBinning","page":"Home","title":"AdaptiveHierarchicalRegularBinning","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for AdaptiveHierarchicalRegularBinning.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [AdaptiveHierarchicalRegularBinning]","category":"page"},{"location":"#AdaptiveHierarchicalRegularBinning.CountSortDetails","page":"Home","title":"AdaptiveHierarchicalRegularBinning.CountSortDetails","text":"struct CountSortDetails{B, D}\n\nContains all of the details needed for the countsort_xxx_impl! algorithms.\n\nFields\n\nB: The length of the bit mask to use.\nD: The leading dimension.\nlo: Lower array bound.\nhi: Higher array bound.\nl: The level of recursion.\n\n\n\n\n\n","category":"type"},{"location":"#AdaptiveHierarchicalRegularBinning.CountSortDetails-Tuple{AdaptiveHierarchicalRegularBinning.RadixSortDetails}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.CountSortDetails","text":"CountSortDetails(rsd)\n\n\nGets a CountSortDetails object from a RadixSortDetails object.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.CountSortDetails-Tuple{Any, Any, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.CountSortDetails","text":"CountSortDetails(bitlen, lo, hi)\n\n\nConstructs a new CountSortDetails object.\n\nArguments\n\nbitlen: The bit length of the bit mask used by the algorithm.\nlo: Lower array bound.\nhi: Higher array bound.\n\nKeyword Arguments\n\ndims: The leading dimension.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.NodeInfo","page":"Home","title":"AdaptiveHierarchicalRegularBinning.NodeInfo","text":"struct NodeInfo\n\nRepresents a node of the tree.\n\nFields\n\nlo: The start index of the represented group.\nhi: The stop index of the represented group.\ndpt: The depth of the node\nnidx: The index of the node.\npidx: The index of the parent node.\n\n\n\n\n\n","category":"type"},{"location":"#AdaptiveHierarchicalRegularBinning.RadixSortDetails","page":"Home","title":"AdaptiveHierarchicalRegularBinning.RadixSortDetails","text":"struct RadixSortDetails{B, D}\n\nContains all of the details needed for the radixsort_xxx_xxx_impl! algorithms.\n\nFields\n\nB: The length of the bit mask to use.\nD: The leading dimension.\ncsd: The underlying CountSortDetails object.\npar_th: Threshold for parallel-to-parallel execution.\nseq_th: Threshold for parallel-to-sequential execution.\nsml_th: Threshold for skiping small arrays.\ndpt_th: recursion threshold.\n\n\n\n\n\n","category":"type"},{"location":"#AdaptiveHierarchicalRegularBinning.RadixSortDetails-Tuple{Any, Any, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.RadixSortDetails","text":"RadixSortDetails(bitlen, lo, hi)\n\n\nConstructs a new RadixSortDetails object.\n\nArguments\n\nbitlen: The bit length of the bit mask used by the algorithm.\nlo: Lower array bound.\nhi: Higher array bound.\n\nKeyword Arguments\n\npar_th: Threshold for parallel-to-parallel execution.\nseq_th: Threshold for parallel-to-sequential execution.\nsml_th: Threshold for skiping small arrays.\ndpt_th: recursion threshold.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.SpatialTree","page":"Home","title":"AdaptiveHierarchicalRegularBinning.SpatialTree","text":"struct SpatialTree{T, E, B, D}\n\nRepresents the tree.\n\nFields\n\nT: The element type of the points.\nE: The element type of the encoded points.\nB: The bit length of each bit group in the encoded space.\nD: The leading dimension.\ninfo: The tree information.\nnidx: Index to the current root of the tree.\n\n\n\n\n\n","category":"type"},{"location":"#AdaptiveHierarchicalRegularBinning.TreeInfo","page":"Home","title":"AdaptiveHierarchicalRegularBinning.TreeInfo","text":"struct TreeInfo{T, E, B, D}\n\nRepresents the tree information.\n\nFields\n\nT: The element type of the points.\nE: The element type of the encoded points.\nB: The bit length of each bit group in the encoded space.\nD: The leading dimension.\npoints: The cloud of points that the tree spans.\nencoded: The encoded cloud of points.\nscale: The original scale of the dimensions of the points.\noffset: The per dimension offset of the points.\nnodes: Per node information.\nchildren: Per node children.\nmaxdepth: The maximum depth of the tree.\nsmlth: Small threshold.\n\n\n\n\n\n","category":"type"},{"location":"#AdaptiveHierarchicalRegularBinning.bit_interleave-Tuple{AbstractVector{<:Unsigned}}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.bit_interleave","text":"bit_interleave(W)\n\nInterleaves an array of unsigned integers.\n\nArguments\n\nW: The array\n\nExamples\n\njulia> bit_interleave([0x00FF, 0x000F]) |> bitstring\n\"0101010111111111\"\n\njulia> bit_interleave([0x000F, 0x00FF]) |> bitstring\n\"1010101011111111\"\n\njulia> bit_interleave([0x0080, 0x0001]) |> bitstring\n\"0100000000000010\"\n\njulia> bit_interleave([0x0001, 0x0080]) |> bitstring\n\"1000000000000001\"\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.bit_space-Tuple{Unsigned, Integer}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.bit_space","text":"bit_space(w, n)\n\nEvenly space out bits.\n\nArguments\n\nw: The word to space out\nn: The amount of padding\n\nExamples\n\njulia> bit_space(0x00ABCDEF, 0) |> bitstring\n\"00000000101010111100110111101111\"\n\njulia> bit_space(0x00ABCDEF, 1) |> bitstring\n\"10100000101000101010100010101010\"\n\njulia> bit_space(0x00ABCDEF, 2) |> bitstring\n\"00000100100100100000100100100100\"\n\njulia> bit_space(0x00ABCDEF, 3) |> bitstring\n\"10001000100000001000100010001000\"\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.bit_space_mask-Tuple{Type{<:Unsigned}, Integer}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.bit_space_mask","text":"bit_space_mask(_, _)\n\n\nGenerates a spaced mask.\n\nArguments\n\nT: The resulting type\nn: The amount of padding\n\nExamples\n\njulia> bit_space_mask(UInt8, 0) |> bitstring\n\"11111111\"\n\njulia> bit_space_mask(UInt8, 1) |> bitstring\n\"01010101\"\n\njulia> bit_space_mask(UInt8, 2) |> bitstring\n\"01001001\"\n\njulia> bit_space_mask(UInt8, 3) |> bitstring\n\"00010001\"\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.bitlen-Union{Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails{B}}, Tuple{B}} where B","page":"Home","title":"AdaptiveHierarchicalRegularBinning.bitlen","text":"Gets the bit length the object.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.bitmask-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.bitmask","text":"bitmask(csd)\n\n\nGets the bit mask the object.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.count_nodes-NTuple{6, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.count_nodes","text":"count_nodes(R, lo, hi, l, depth, bitlen)\n\n\nCounts all nodes in the tree.\n\nArguments\n\nR: The array to search.\nlo: The start index.\nhi: The stop index.\nl: Maximum depth.\ndepth: The current depth in the tree. (Starts at 0)\nbitlen: The bitlen of the each bit group.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.countsort_par_impl!-Union{Tuple{TC}, Tuple{TI}, Tuple{TR}, Tuple{TV}, Tuple{TV, TR, TI, TV, TR, TI, TC, AdaptiveHierarchicalRegularBinning.CountSortDetails}} where {TV<:AbstractArray, TR<:(AbstractVector{<:Unsigned}), TI<:(AbstractVector{<:Unsigned}), TC<:(AbstractMatrix{<:Unsigned})}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.countsort_par_impl!","text":"countsort_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, csd)\n\n\nParallel countsort.\n\nRemarks\n\nInputs Va, Ra, Ia are not altered.\n\nArguments\n\nVa: The cloud of points.\nRa: The morton transformation for each point of Va\nIa: The current permutation of Va.\nVb: A matrix to store the sorted version of Va.\nRb: A vector to store the sorted version of Ra.\nIb: A vector to store the resulting permutation.\nC: A matrix acting as a working space for the countsort algorithm.\ncsd: The CountSortDetails object used to configure the algorithm.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.countsort_seq_impl!-Union{Tuple{TC}, Tuple{TI}, Tuple{TR}, Tuple{TV}, Tuple{TV, TR, TI, TV, TR, TI, TC, AdaptiveHierarchicalRegularBinning.CountSortDetails}} where {TV<:AbstractArray, TR<:(AbstractVector{<:Unsigned}), TI<:(AbstractVector{<:Unsigned}), TC<:(AbstractVector{<:Unsigned})}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.countsort_seq_impl!","text":"countsort_seq_impl!(Va, Ra, Ia, Vb, Rb, Ib, C, csd)\n\n\nSequential countsort.\n\nRemarks\n\nInputs Va, Ra, Ia are not altered.\n\nArguments\n\nVa: The cloud of points.\nRa: The morton transformation for each point of Va\nIa: The current permutation of Va.\nVb: A matrix to store the sorted version of Va.\nRb: A vector to store the sorted version of Ra.\nIb: A vector to store the resulting permutation.\nC: A vector acting as a working space for the countsort algorithm.\ncsd: The CountSortDetails object used to configure the algorithm.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.depth-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.depth","text":"depth(csd)\n\n\nGets the level of recursion of the object.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.get_node_hi-NTuple{5, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.get_node_hi","text":"get_node_hi(R, lo, hi, depth, bitlen)\n\n\nLinear searches the next node high.\n\nArguments\n\nR: The array to search.\nlo: The start index.\nhi: The stop index.\ndepth: The current depth in the tree. (Starts at 0)\nbitlen: The bitlen of the each bit group.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.high-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.high","text":"high(csd)\n\n\nGets the higher index of the object.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.leaddim-Union{Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails{<:Any, D}}, Tuple{D}} where D","page":"Home","title":"AdaptiveHierarchicalRegularBinning.leaddim","text":"leaddim(_)\n\n\nGets the leading dimension the object.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.low-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.low","text":"low(csd)\n\n\nGets the lower index of the object.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.make_tree-NTuple{7, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.make_tree","text":"make_tree(V, R, maxdpt, smlth, bitlen, scale, offset)\n\n\nCreates a tree representation of a Morton Array.\n\nArguments\n\nV: Cloud of points.\nR: The array to convert.\nmaxdpt: Maximum depth.\nsmlth: Small threshold.\nscale: Scalar for the original coordinates.\noffset: Per dimension offset.\n\nKeyword Arguments\n\ndims: Leading dimension\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.make_tree_impl-Tuple{Any, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.make_tree_impl","text":"make_tree_impl(node, idx)\n\n\nCreates a tree representation of a Morton Array.\n\nArguments\n\nnode: The current node to build.\nidx: Current node index.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.next-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails, Any, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.next","text":"next(csd, lo, hi)\n\n\nConstructs a new CountSortDetails object with a deeper level of recursion.\n\nArguments\n\nlo: New lower array bound.\nhi: New higher array bound.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.next-Tuple{AdaptiveHierarchicalRegularBinning.RadixSortDetails, Any, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.next","text":"next(rsd, lo, hi)\n\n\nConstructs a new RadixSortDetails object with a deeper level of recursion.\n\nArguments\n\nlo: New lower array bound.\nhi: New higher array bound.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.quantize-Tuple{Any, Any, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.quantize","text":"quantize(T, x, l)\n\n\nQuantizes a single coordinate.\n\nArguments\n\nT: The resulting type.\nx: The input coordinate.\nl: The levels of the quantizer.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.radixsel-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.radixsel","text":"radixsel(csd, x)\n\n\nSelects the corresponding radix.\n\nArguments\n\ncsd: The CountSortDetails object to base the selection on.\nx: The value from which the radix is computed.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.radixshft-Tuple{Any, Any, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.radixshft","text":"radixshft(x, depth, bitlen)\n\n\nShifts the corresponding amount of bits.\n\nArguments\n\nx: The value from which the radix is computed.\ndepth: The current depth.\nbitlen: The bit length of the bit groups.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.radixsort_par_par_impl!-Union{Tuple{TI}, Tuple{TA}, Tuple{TR}, Tuple{TV}, Tuple{TV, TR, TI, TV, TR, TI, AbstractVector{Bool}, AdaptiveHierarchicalRegularBinning.RadixSortDetails, AdaptiveHierarchicalRegularBinning.Allocator{TA}}} where {TV<:(AbstractMatrix), TR<:(AbstractVector{<:Unsigned}), TA<:Unsigned, TI<:AbstractVector{TA}}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.radixsort_par_par_impl!","text":"radixsort_par_par_impl!(Va, Ra, Ia, Vb, Rb, Ib, P, rsd, allocator)\n\n\nParallel-to-Parallel radixsort.\n\nArguments\n\nVa: The cloud of points.\nRa: The morton transformation for each point of Va\nIa: The current permutation of Va.\nVb: An auxiliary matrix to store the partially sorted version of Va.\nRb: An auxiliary vector to store the partially sorted version of Ra.\nIb: An auxiliary vector to store the resulting permutation.\nP: A ping-pong matrix that denotes where the latest results are.\nrsd: The RadixSortDetails object used to configure the algorithm.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.radixsort_par_seq_impl!-Union{Tuple{TI}, Tuple{TA}, Tuple{TR}, Tuple{TV}, Tuple{TV, TR, TI, TV, TR, TI, AbstractVector{Bool}, AdaptiveHierarchicalRegularBinning.RadixSortDetails, AdaptiveHierarchicalRegularBinning.Allocator{TA}}} where {TV<:(AbstractMatrix), TR<:(AbstractVector{<:Unsigned}), TA<:Unsigned, TI<:AbstractVector{TA}}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.radixsort_par_seq_impl!","text":"Parallel-to-Sequential radixsort.\n\nArguments\n\nVa: The cloud of points.\nRa: The morton transformation for each point of Va\nIa: The current permutation of Va.\nVb: An auxiliary matrix to store the partially sorted version of Va.\nRb: An auxiliary vector to store the partially sorted version of Ra.\nIb: An auxiliary vector to store the resulting permutation.\nP: A ping-pong matrix that denotes where the latest results are.\nrsd: The RadixSortDetails object used to configure the algorithm.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.radixsort_seq_seq_impl!-Union{Tuple{TI}, Tuple{TA}, Tuple{TR}, Tuple{TV}, Tuple{TV, TR, TI, TV, TR, TI, AbstractVector{Bool}, AdaptiveHierarchicalRegularBinning.RadixSortDetails, AdaptiveHierarchicalRegularBinning.Allocator{TA}}} where {TV<:(AbstractMatrix), TR<:(AbstractVector{<:Unsigned}), TA<:Unsigned, TI<:AbstractVector{TA}}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.radixsort_seq_seq_impl!","text":"Sequential-to-Sequential radixsort.\n\nArguments\n\nVa: The cloud of points.\nRa: The morton transformation for each point of Va\nIa: The current permutation of Va.\nVb: An auxiliary matrix to store the partially sorted version of Va.\nRb: An auxiliary vector to store the partially sorted version of Ra.\nIb: An auxiliary vector to store the resulting permutation.\nP: A ping-pong matrix that denotes where the latest results are.\nrsd: The RadixSortDetails object used to configure the algorithm.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.regural_bin-NTuple{4, Any}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.regural_bin","text":"regural_bin(RT, V, l, smlth)\n\n\nConstructs the tree.\n\nArguments\n\nRT: The type of the morton vector.\nV: The cloud of points.\nl: The maximum tree depth.\nsmlth: Small threshold.\n\nKeyword Arguments\n\ndims: Leading dimension\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.spatial_encode!-Union{Tuple{D}, Tuple{AbstractVector{<:Integer}, AbstractMatrix, Integer}} where D","page":"Home","title":"AdaptiveHierarchicalRegularBinning.spatial_encode!","text":"spatial_encode!(R, V, l)\n\n\nEncodes a set of points using the morton encoding.\n\nArguments\n\nR: The resulting vector.\nV: The input matrix. Will be scaled and translated.\nl: The levels of the encoding.\n\nKeyword Arguments\n\ndims: The leading dimension.\ncenter: Whether to divide the slack evenly.\n\n\n\n\n\n","category":"method"},{"location":"#AdaptiveHierarchicalRegularBinning.translate_scale_vals-Tuple{AbstractMatrix}","page":"Home","title":"AdaptiveHierarchicalRegularBinning.translate_scale_vals","text":"translate_scale_vals(V)\n\n\nCompute displacement vector and the scale to transform the cloud to a unit hypercube.\n\nArguments\n\nV: The input matrix.\n\nKeyword Arguments\n\ndims: The leading dimension.\ncenter: Whether to divide the slack evenly.\n\n\n\n\n\n","category":"method"},{"location":"#Base.eltype-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails}","page":"Home","title":"Base.eltype","text":"eltype(csd)\n\n\nThe index type of the struct.\n\n\n\n\n\n","category":"method"},{"location":"#Base.length-Tuple{AdaptiveHierarchicalRegularBinning.CountSortDetails}","page":"Home","title":"Base.length","text":"length(csd)\n\n\nGets the length of the object.\n\n\n\n\n\n","category":"method"}]
}
