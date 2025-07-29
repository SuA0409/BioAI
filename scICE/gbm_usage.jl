run(`git clone https://github.com/Mathbiomed/scICE`)
cd("scICE")

import Pkg
Pkg.activate(".")
Pkg.instantiate()

ENV["NUM_CORES"] = "12"

using CUDA, CSV, DataFrames, scLENS
include("C:/Users/user/scICE/src/scICE.jl")
using CairoMakie
CairoMakie.activate!(type="png")


cur_dev = if CUDA.has_cuda()
    "gpu"
else
    "cpu"
end

# preprocessing input data
using Pkg
Pkg.add(["CodecZlib", "CSV", "DataFrames", "MatrixMarket", "SparseArrays"])
using CodecZlib, MatrixMarket, SparseArrays

path = raw"C:/Users/user/Desktop/gbm_data/GBM2_extracted/GBM2_filtered_feature_bc_matrix"

# barcodes.tsv.gz 읽기
barcodes_file = joinpath(path, "barcodes.tsv.gz")
open(GzipDecompressorStream, barcodes_file) do io
    barcodes = readlines(io)
    global barcodes_vec = barcodes
end

# features.tsv.gz 읽기 (두 번째 컬럼만 추출: gene 이름)
features_file = joinpath(path, "features.tsv.gz")
open(GzipDecompressorStream, features_file) do io
    features_df = CSV.read(io, DataFrame; header=false, delim='\t')
    global gene_names = features_df[:, 2]  # 두 번째 컬럼
end

# matrix.mtx.gz 읽기 (CSR 희소 행렬로 변환)
matrix_file = joinpath(path, "matrix.mtx.gz")
mtx = MatrixMarket.mmread(matrix_file)
mtx = convert(SparseMatrixCSC{Float64, Int}, mtx)

println("Original matrix shape: ", size(mtx))  # (genes, cells)

# 보통 (genes x cells) 형태 → transpose 필요
mtx_t = transpose(mtx)
println("Transposed matrix shape: ", size(mtx_t))  # (cells, genes)

# DataFrame으로 변환 (주의: 매우 클 수 있음!)
pre_df = DataFrame(mtx_t, :auto)
rename!(pre_df, Symbol.(gene_names); makeunique=true)
# features.tsv.gz 안에 같은 이름의 유전자가 여러 개 존재 -> TBCE, TBCE_1, TBCE_2 같은 식으로 중복 이름에 숫자 붙여 고유 컬럼명으로 만듦.
pre_df[!, :cell] = barcodes_vec
select!(pre_df, :cell, Not(:cell))  # 'cell'을 첫 컬럼으로 이동

println(size(pre_df))
display(pre_df)

# PCA
sclens_embedding = scLENS.sclens(pre_df, device_=cur_dev)
CSV.write("out/pca.csv", sclens_embedding[:pca_n1])

scLENS.apply_umap!(sclens_embedding)
CSV.write("out/umap.csv", DataFrame(sclens_embedding[:umap], :auto))

panel_0 = scLENS.plot_embedding(sclens_embedding)
save("out/umap_dist.png",panel_0)

clustering!(sclens_embedding)
# By default, scICE explores cluster numbers ranging from 1 to 20 
# (this is the default value for the optional second argument r, 
# as seen in the function signature clustering!(a_dict, r=[1,20]; ...)). 
# If you wish to focus the analysis on a specific range of cluster numbers, 
# for instance, from 5 to 10 clusters, you provide this range as the second argument:
clustering!(sclens_embedding, [5,10])
# This enables you to find consistent cluster labels more efficiently within an anticipated range.

# 9. Inconsistency Coefficient Visualization: Visualize the Inconsistency coefficient and save it
panel_1 = plot_ic(sclens_embedding)
save("out/ic_plot.png",panel_1)

# 10. Consistent Cluster Label Extraction: Extract consistent cluster labels using get_rlabel! and save them to a CSV file. 
# This function filters labels based on an Inconsistency Coefficient (IC) threshold.
label_out = get_rlabel!(sclens_embedding)
CSV.write("out/consistent_labels.csv", label_out)
# The IC threshold parameter (th) defaults to 1.005. This value is passed as the 
# optional second argument to get_rlabel! and can be adjusted if needed. For example, 
# to change the threshold to 1.01, you would call the function like this:
label_out = get_rlabel!(sclens_embedding, 1.01)

# 11. Cluster Visualization: Set the number of clusters and visualize them with labels
n_clusters = 8
panel_2 = scLENS.plot_embedding(sclens_embedding, label_out[!, "l_$n_clusters"])
save("out/umap_dist_with_label$n_clusters.png",panel_2)

# 12. Save result as AnnData
scLENS.save_anndata("out/test.h5ad",sclens_embedding)
