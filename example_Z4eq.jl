# Clone the repository
run(`git clone https://github.com/Mathbiomed/scICE`)
cd("scICE")

# Active the project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

# 1. Configure Processing Cores: Set the number of CPU cores for processing
ENV["NUM_CORES"] = "12"

# 2. Set up the environment: Load the necessary packages and include the local scICE file
using CUDA, CSV, DataFrames, scLENS
# pwd() 실행 시 -> "C:\\Users\\user\\scICE"
include("src/scICE.jl")
using CairoMakie
CairoMakie.activate!(type="png")

# 3. Device Selection: The device (CPU or GPU) is automatically selected based on CUDA availability
cur_dev = if CUDA.has_cuda()
    "gpu"
else
    "cpu"
end

# 4. Data Preprocessing: Load your single-cell data (example uses compressed CSV) and preprocess it
ndf = scLENS.read_file(raw"C:/Users/user/Desktop/Z4eq.csv.gz")
pre_df = scLENS.preprocess(ndf)

# 5. Embedding Creation: Create an embedding for the preprocessed data using scLENS
sclens_embedding = scLENS.sclens(pre_df, device_=cur_dev)
CSV.write("out/pca.csv", sclens_embedding[:pca_n1])

# 6. UMAP Transformation: Apply UMAP to the embedding and save the results
scLENS.apply_umap!(sclens_embedding)
CSV.write("out/umap.csv", DataFrame(sclens_embedding[:umap], :auto))

# 7. Visualization: Plot the UMAP distribution and save the output
panel_0 = scLENS.plot_embedding(sclens_embedding)
save("out/umap_dist.png",panel_0)

# 8. Applying scICE: Apply scICE clustering to the embedding
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