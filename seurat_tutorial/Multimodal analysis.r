####### Multimodal analysis #######

# 1. Load in the data
library(Seurat)
library(ggplot2)
library(patchwork)
# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative controls for the protein measurements. 
# For this reason, the gene expression matrix has HUMAN_ or MOUSE_ appended to the beginning of each gene.
cbmc.rna <- as.sparse(read.csv(file = "C:/Users/user/Desktop/Seurat/tutorial/2/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
    sep = ",", header = TRUE, row.names = 1))

# To make life a bit easier going forward, we're going to discard all but the top 100 most highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix.
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "C:/Users/user/Desktop/Seurat/tutorial/2/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
    sep = ",", header = TRUE, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical column names.
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))


# 2. Setup a Seurat object, add the RNA and protein data
# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)
# We can see that by default, the cbmc object contains an assay storing RNA measurement
Assays(cbmc)

# create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = cbmc.adt)
# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay
# Validate that the object now contains multiple assays
Assays(cbmc)

# Extract a list of features measured in the ADT assay
rownames(cbmc[["ADT"]])

# Note that we can easily switch back and forth between the two assays to specify the default for visualization and analysis.
# List the current default assay
DefaultAssay(cbmc)

# Switch the default to ADT
DefaultAssay(cbmc) <- "ADT"
DefaultAssay(cbmc)


# 3. Cluster cells on the basis of their scRNA-seq profiles
# quick clustering of the PBMCs based on the scRNA-seq data
# Note that all operations below are performed on the RNA assay Set and verify that the default assay is RNA.
DefaultAssay(cbmc) <- "RNA"
DefaultAssay(cbmc)

# Guide tutorial -2,700 PBMCs 참고
# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)   # Normalizing the data (line 51)
cbmc <- FindVariableFeatures(cbmc)   # Identification of highly variable features(feature selection) (line 57)
cbmc <- ScaleData(cbmc)   # Scaling the data (line 67)
cbmc <- RunPCA(cbmc, verbose = FALSE)   # Linear dimensional reduction(e.g. PCA) (line 77)
cbmc <- FindNeighbors(cbmc, dims = 1:30)   # Cluster the cells (line 96)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)   # Cluster the cells; 해상도 매개변수 포함 (line 99)
cbmc <- RunUMAP(cbmc, dims = 1:30)   # Non-linear dimensional reduction(UMAP/tSNE) (line 105)
DimPlot(cbmc, label = TRUE)   # individual clusters (line 109)


# 4. Visualize multiple modalities side-by-side
# Normalize ADT data,
DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2)
DefaultAssay(cbmc) <- "RNA"

# Note that the following command is an alternative but returns the same result
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")

# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can visualize one or the other.
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")
# place plots side-by-side
p1 | p2

# Alternately, we can use specific assay keys to specify a specific modality Identify the key for the RNA and protein assays.
Key(cbmc[["RNA"]])
Key(cbmc[["ADT"]])
# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(cbmc, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2
    # Key를 사용하는 이유: 
    # 여러 modality(e.g. RNA, ADT protein 등)의 feature 이름(e.g. CD19)이 RNA와 ADT에 동시에 존재할 수 있음.
    # Key를 사용하여 서로 다른 이름으로 관리되도록 자동으로 prefix를 붙여줌. 
    # => key를 사용하면 어떤 assay에서 온 feature인지 명확하게 구분 가능.


# 5. Identify cell surface markers for scRNA-seq clusters
# as we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the surface
VlnPlot(cbmc, "adt_CD19")

# we can also identify alternative protein and RNA markers for this cluster through differential expression
adt_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "ADT")
rna_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "RNA")
head(adt_markers)
head(rna_markers)


# 6. Additional visualizations of multimodal data
# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells
# if desired by using HoverLocator and FeatureLocator
FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
# view relationship between protein and RNA
FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")
# Let's look at the raw (non-normalized) ADT counts. You can see the values are quite high,
# particularly in comparison to RNA values. This is due to the significantly higher protein copy number in cells, 
# which significantly reduces 'drop-out' in ADT data
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")


# + Loading data from 10X multi-modal experiments
### 데이터 다운로드를 못해서 실행 X ###
pbmc10k.data <- Read10X(data.dir = "/brahms/shared/vignette-data/pbmc10k/filtered_feature_bc_matrix/")
# 전처리: 단백질 항목 이름 정리
rownames(x = pbmc10k.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
    x = rownames(x = pbmc10k.data[["Antibody Capture"]]))
# e.g. "CD19_control_TotalSeqB" → "CD19"
#       "CD3_TotalSeqB" → "CD3"
# 이후 분석이나 시각화에서 feature1 = "adt_CD19" 등 사용시 이름이 정확히 매칭되어 오류 안 남.

pbmc10k <- CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 3, min.features = 200)
pbmc10k <- NormalizeData(pbmc10k)
pbmc10k[["ADT"]] <- CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])
pbmc10k <- NormalizeData(pbmc10k, assay = "ADT", normalization.method = "CLR")

plot1 <- FeatureScatter(pbmc10k, feature1 = "adt_CD19", feature2 = "adt_CD3", pt.size = 1)
plot2 <- FeatureScatter(pbmc10k, feature1 = "adt_CD4", feature2 = "adt_CD8a", pt.size = 1)
plot3 <- FeatureScatter(pbmc10k, feature1 = "adt_CD3", feature2 = "CD3E", pt.size = 1)
(plot1 + plot2 + plot3) & NoLegend()

plot <- FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3") + NoLegend() + theme(axis.title = element_text(size = 18),
    legend.text = element_text(size = 18))
ggsave(filename = "../output/images/citeseq_plot.jpg", height = 7, width = 12, plot = plot, quality = 50)
