# Analyze the reference transcriptome of 10 bead-enriched subpopulations of
# PBMCs from Donor A. See Zheng, 2017 fore more details.

using PythonCall, Tar, CodecZlib, LinearAlgebra, SparseArrays

const scanpy = pyimport("scanpy")
const scipy = pyimport("scipy")
const numpy = pyimport("numpy")
const pandas = pyimport("pandas")

function download10xgenomics(url::String)::String
  tarball = download( url )
  tar_gz = open( tarball )
  tar = GzipDecompressorStream( tar_gz )
  dir = Tar.extract( tar )
  close( tar )

  return dir
end

"""
    pysparseconvert( X::Py )::SparseMatrixCSC

Convert a sparse matrix from Python to Julia.
"""
function pysparseconvert( X::Py )::SparseMatrixCSC
  (ii,jj,vv) = scipy.sparse.find(X)
  ii = pyconvert( Vector, ii ) .+ 1
  jj = pyconvert( Vector, jj ) .+ 1
  vv = pyconvert( Vector, vv )
  (m, n) = pyconvert( Tuple, X.shape )

  return sparse( ii, jj, vv, m, n )
end

##################################################
# # Prepare scanpy
# figdir = mkpath( joinpath( dirname( @__FILE__ ), "figures_zheng_references" ) )
scanpy.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
scanpy.logging.print_header()
# scanpy.settings.figdir = figdir
# scanpy.settings.set_figure_params(dpi=80, facecolor="white")
scanpy.settings.autosave = false
scanpy.settings.autoshow = false

# # Download datasets
# We gather each reference transcriptome from the 10x genomics repository and
# combine them into a single matrix.
urls = Dict(
    "CD14+ Monocytes" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/cd14_monocytes_filtered_gene_bc_matrices.tar.gz",
    "CD19+ B Cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz",
    "CD4+ Helper T Cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd4_t_helper/cd4_t_helper_filtered_gene_bc_matrices.tar.gz",
    "CD4+/CD25+ Regulatory T Cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/regulatory_t/regulatory_t_filtered_gene_bc_matrices.tar.gz",
    "CD4+/CD45RA+/CD25- Naive T cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_t/naive_t_filtered_gene_bc_matrices.tar.gz",
    "CD4+/CD45RO+ Memory T Cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/memory_t/memory_t_filtered_gene_bc_matrices.tar.gz",
    "CD56+ Natural Killer Cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd56_nk/cd56_nk_filtered_gene_bc_matrices.tar.gz",
    "CD8+ Cytotoxic T cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cytotoxic_t/cytotoxic_t_filtered_gene_bc_matrices.tar.gz",
    "CD8+/CD45RA+ Naive Cytotoxic T Cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_cytotoxic/naive_cytotoxic_filtered_gene_bc_matrices.tar.gz",
    "CD34+ Cells" => "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd34/cd34_filtered_gene_bc_matrices.tar.gz"  # 45% purity (maybe we shouldn't include it?)
)

strCellTypes = collect(keys(urls));   # the cell types

## Download the transcriptomes for each cell type
adata = map( strCellTypes ) do cell
    dir = download10xgenomics( urls[cell] )
    scanpy.read_10x_mtx(
        joinpath( dir, "filtered_matrices_mex/hg19/" ),  # the directory with the `.mtx` file
        var_names="gene_symbols",                        # use gene symbols for the variable names (variables-axis index)
        cache=false)

end

# after concatenation, the refernce configuration is the "batch" parameter
conc = adata[1].concatenate( adata[2:end]... ) # concatenate the results
conc.var_names_make_unique()  # make sure the variables (feature dimensions) are unique

# # Preprocessing
# The preprocessing steps here follow the normalization and filtering as of
# Zheng, 2017.

n_top_genes = 1000                            # how many highly-variable genes to keep
scanpy.pp.filter_genes(conc, min_counts=1)    # only consider genes with more than 1 count
scanpy.pp.normalize_per_cell(                 # normalize with total UMI count per cell
     conc, key_n_counts="n_counts_all"
)
filter_result = scanpy.pp.filter_genes_dispersion(  # select highly-variable genes
    conc.X, flavor="cell_ranger", n_top_genes=n_top_genes, log=false
)

# subset the genes
conc = @pyeval (conc = conc, filter_result = filter_result) => `conc[:, filter_result.gene_subset]`

scanpy.pp.normalize_per_cell(conc)            # renormalize after filtering
scanpy.pp.log1p(conc)                         # log transform: conc.X = log(conc.X + 1)
scanpy.pp.scale(conc)                         # scale to unit variance and shift to zero mean

# ## Principal component analysis
# Reduce the dimensionality of the data by running principal component analysis
# (PCA), which reveals the main axes of variation and denoises the data.

scanpy.tl.pca(conc, svd_solver="arpack", n_comps = 50);

########################################################################################
## Experiments:
# Consider as point coordinates the PCA of the above data

X = pyconvert(Matrix,conc.obsm["X_pca"])

# 94655×50 Matrix{Float32}
# The dimensions are too many, but are ordered in descending order
# > sqrt.(sum(X.^2,dims=1))
#
# Use the first 30 columns as point coordinates. This should allow us to get a Level == 4
# tree with 128 bit integer encodings.
# Questions:
#  - What is the distribution of points per leaf box for B == 2000?
#  - What percentage of the NxN distances are needed to find k-NN for k = 2.^3:6?
#  - What if I repeat the above for X[:,1:50] and two levels and X[:,1:40] and 3 levels?

using AdaptiveHierarchicalRegularBinning

# TODO: Accept X as the input but pass 1:30 as the point-coords-idx param
cols = 1:30
PXT = X[:, cols] |> transpose |> collect
dpt = 4
smlth = 2000
t = AdaptiveHierarchicalRegularBinning.regural_bin(UInt128, PXT, dpt, smlth; dims=2);

include("knn.jl")
k = 2^3
indices, distances, levels = knn(t, PXT, k);
