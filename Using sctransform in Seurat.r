####### Using sctransform in Seurat #######
library(Seurat)
library(ggplot2)
library(sctransform)


# Load data and create Seurat object
pbmc_data <- Read10X(data.dir = "C:/Users/user/Desktop/Seurat/tutorial/1/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)


# Apply sctransform normalization
# store mitochondrial percentage in object meta data
# "^MT-"로 시작하는 유전자의 발현 비율 계산해서 "percent.mt"라는 메타데이터 컬럼에 저장
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
# "percet.mt"를 회귀변수로 넣어 정규화하고 변환(mito effect 제거)
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


# Perform dimensionality reduction by PCA and UMAP embedding
# These are now standard steps in the Seurat workflow for visualization and clustering
# SCTransform으로 정규하된 데이터를 기반으로,
# PCA -> UMAP -> 이웃 찾기 -> 클러스터링 -> 시각화
# PCA(주성분 분석)으로 각 셀을 고차원에서 저차원 공간으로 축소해 유사한 셀들이 가까이 위치하도록 만듦
pbmc <- RunPCA(pbmc, verbose = FALSE)
# PCA 결과의 1~30번째 주성분을 이용해 UMAP 수행 (셀 간 유사도를 바탕으로 2D 시각화)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
# 각 셀이 다른 셀들과 얼마나 가까운지 이웃 관계 계산 (KNN 그래프)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# 이웃 정보 그래프를 기반으로 그래프 클러스터링 수행.
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE)

# These are now standard steps in the Seurat workflow for visualization and clustering
# Visualize canonical marker genes as violin plots.
VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), pt.size = 0.2, ncol = 4)
# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2, ncol = 3)
FeaturePlot(pbmc, features = c("CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"), pt.size = 0.2, ncol = 3)
