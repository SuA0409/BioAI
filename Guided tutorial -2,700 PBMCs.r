install.packages('Seurat')
library(Seurat)

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)



####### Guided tutorial -2,700 PBMCs #######
# 10X Genomics에서 제공하는 PBMC 데이터 세트 분석

library(dplyr)
library(Seurat)
library(patchwork)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/user/Desktop/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
    # CreateSeuratObject()에서 고유 유전자 수(unique genes)와 총 분자수(total molecules) 자동으로 계산됨.
pbmc
pbmc.data[c("CD3D","TCL1A","MS4A1"), 1:30]   # 카운트 행렬의 데이터 출력값 형식

# Preprocessing (전처리)
 # QC and selecting cells for further analysis (QC 및 추가 분석을 위한 세포 선택)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    # 미토콘드리아 QC 메트릭 계산
    # calculate mitochondrial QC metrics with the PercentageFeatureSet(),
        # which calculates the percentage of counts originating from a set of features
    # MT- : 미토콘드리아 유전자 세트로 시작하는 모든 유전자 세트 사용
head(pbmc@meta.data, 5)   # Seurat의 QC 지표 출력
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# 2500개 이상 또는 200개 미만의 고유 기능 수를 갖는 셀 + 미토콘드리아 수가 5% 이상인 세포 필터링
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# 원하는 데이터만 추출

 # Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    # LogNormalize : 각 셀의 특징 표현 측정값을 전체 표현값으로 정규화 -> 스케일 인자(기본 10,000)을 곱함 -> 로그 변환
    # pbmc <- NormalizeData(pbmc) : 이렇게 해도 동일함 (기본값으로 제공하기 때문)
    # 사용 전제조건: 각 세포가 원래 동일한 수의 RNA 분자 포함

 # Identification of highly variable features (feature selection)
 # 데이터 세트에서 어떤 세포에서는 발현량이 높고 다른 세포에서는 발현량이 낮은 특징들 계산
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

 # Scaling the data
 # 선형변환(스케일링) -> 차원축소기법(e.g. PCA)
 # Normalization을 진행해야만 본 단계가 진행됨.
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
    # 세포 전체의 평균 발현 0, 세포 간 분산 1이 되도록 함. (동일한 가중치를 부여하므로 고도로 발현된 유전자가 지배적 X)
    # default: variable features만 스케일됨. -> features 인자에서 추가해서 스케일하면 됨.
 # 원치 않는 변동 원인 제거방법 (e.g. 세포 주기 단계나 미토콘드리아 오염)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

 # Linear dimensional reduction (e.g. PCA)
 # default: 이전에 결정된 변수 특성만 입력으로 사용됨. -> ScaleData()의 features 인자에서 정의한 subset 사용.
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    # 가장 많은 양(+)의 로딩과 음(-)의 로딩을 갖는 유전자 목록 출력
    # 데이터 세트의 단일 세포에서 상관관계(또는 역상관관계)를 보이는 유전자 모듈
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)   # 그 중 상위 5개만 출력
# Visualizing both cells and features that define the PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
    # DimHeatmap(): 이질성의 주요 원인 쉽게 탐색 가능. 후속 분석에 포함할 PC를 결정할 때 유용.
    # PCA 점수에 따라 cell과 feature가 정렬됨.

# Determine the 'dimensionality' of the dataset
ElbowPlot(pbmc)
    # 각 주성분이 설명하는 분산의 백분율을 기준으로 주성분의 순위를 매김.

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
    # PCA 공간에서 유클리드 거리를 기반으로 KNN 그래프 구성
    # 두 셀의 로컬 이웃(자카르 유사도)에서 공통적으로 겹치는 부분을 기반으로 두 셀 간의 에지 가중치를 세분화
pbmc <- FindClusters(pbmc, resolution = 0.5)
    # 해상도 매개변수 포함. (값이 클수록 클러스터 수 증가)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
    # Idents()로 클러스터 찾을 수 있음 (데이터셋이 클수록 최적 해상도 높아지는 경우 많음)

# Non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")   # 객체 저장 가능 (다시 로드해서 시각화)

# Finding differentially expressed(DE) features (cluster biomarker)
# 차등발현(DE)을 통해 클러스터를 정의하는 마커 발굴 
# (ident.1으로 명시된) 단일 클러스터의 양성 및 음성 마커를 다른 모든 세포와 비교하여 식별
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
    # 매개변수를 사용하여 설정할 수 있는 여러 가지 차등 발현 검정
    # e.g. ROC 검정은 개별 마커에 대한 '분류력'을 반환 (0: 무작위, 1: 완벽)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
    # 클러스터 전체의 발현 확률 분포를 보여줌
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
    # tSNE 또는 PCA 플롯에서 특징 발현을 시각화
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
    # 주어진 셀과 특징에 대한 표현형 히트맵 생성
    # e.g. 이 경우, 각 클러스터의 상위 20개 마커 표시 (20개 미만인 경우 모든 마커 표시)

# Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")   # 저장 가능