run(`git clone https://github.com/Mathbiomed/scICE`)
cd("scICE")

import Pkg
Pkg.activate(".")
Pkg.instantiate()

ENV["NUM_CORES"] = "12"

using CUDA, scLENS
include("C:/Users/user/scICE/src/scICE.jl")
using CairoMakie
CairoMakie.activate!(type="png")

cur_dev = if CUDA.has_cuda()
    "gpu"
else
    "cpu"
end

using CSV
using DataFrames
pre_df = CSV.read("C:/Users/user/Desktop/T1D/beta/t2d_beta_treatment_preprocess.csv", DataFrame)

# PCA
sclens_embedding = scLENS.sclens(pre_df, device_=cur_dev)
CSV.write("C:/Users/user/Desktop/T1D/beta/out/pca.csv", sclens_embedding[:pca_n1])

scLENS.apply_umap!(sclens_embedding)
CSV.write("C:/Users/user/Desktop/T1D/beta/out/umap.csv", DataFrame(sclens_embedding[:umap], :auto))

panel_0 = scLENS.plot_embedding(sclens_embedding)
save("C:/Users/user/Desktop/T1D/beta/out/umap_dist.png",panel_0)

clustering!(sclens_embedding, [3, 13])

# Inconsistency Coefficient Visualization: Visualize the Inconsistency coefficient and save it
panel_1 = plot_ic(sclens_embedding)
save("C:/Users/user/Desktop/T1D/beta/out/ic_plot.png",panel_1)

# Consistent Cluster Label Extraction: Extract consistent cluster labels using get_rlabel! and save them to a CSV file. 
# This function filters labels based on an Inconsistency Coefficient (IC) threshold.
label_out = get_rlabel!(sclens_embedding)
CSV.write("C:/Users/user/Desktop/T1D/beta/out/consistent_labels.csv", label_out)
# # The IC threshold parameter (th) defaults to 1.005. This value is passed as the 
# # optional second argument to get_rlabel! and can be adjusted if needed. For example, 
# # to change the threshold to 1.01, you would call the function like this:
# label_out = get_rlabel!(sclens_embedding, 1.01)

# Cluster Visualization: Set the number of clusters and visualize them with labels
n_clusters = 4
panel_2 = scLENS.plot_embedding(sclens_embedding, label_out[!, "l_$n_clusters"])
save("C:/Users/user/Desktop/T1D/beta/out/umap_dist_with_label$n_clusters.png",panel_2)

# Save result as AnnData
scLENS.save_anndata("C:/Users/user/Desktop/T1D/beta/out/test.h5ad", sclens_embedding; device_="cpu")  # 저장 안됨.