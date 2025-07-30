####### Introduction to scRNA-seq integration #######

library(Seurat)
library(SeuratData)
library(patchwork)
# install dataset
install.packages("C:/Users/user/Desktop/Seurat/tutorial/5/ifnb.SeuratData_3.1.0.tar.gz", repos = NULL, type = "source")
# 인터페론 자극세포와 대조군의 두 가지 조건에서 얻은 인간 PBMC 데이터
# 목표: 세포 하위 집단을 공동으로 식별하고, 각 그룹이 조건에 따라 어떻게 다른가.
# load dataset
ifnb <- LoadData("ifnb")
# split the RNA measurements into two layers one for control cells, one for stimulated cells
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb


# 통합 없이 분석만 수행
# Perform analysis without integration
# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# seurat_clusters : 숫자 클러스터 나뉨. (어떤 클러스터에는 STIM이, 어떤 클러스터에는 CTRL이 많이 포함됨.)
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))
# => 조건이 다른 두 샘플(Ctrl vs. STIM)은 같은 세포라도 조건 때문에 발현 패턴이 달라지기 때문에
#   "같은 세포 타입이 조건이 달라도 같은 위치에 모이도록 하여 조건 간 차이를 비교"


# 통합 수행
# Perform integration
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
# Seurat_annotatioons : 진짜 세포 cell type으로 라벨
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
# 같은 세포끼리 잘 뭉쳐 있으므로, 클러스터 자체가 조건과 무관하게 잘 정의되어 조건 간 비교 가능
# => 자극 효과를 세포 타입 별로 비교
DimPlot(ifnb, reduction = "umap", split.by = "stim")


# 보존된 세포 유형 마커 식별
# Identify conserved cell type markers
Idents(ifnb) <- "seurat_annotations"
# FindConservedMarkers() : 동일한 세포 타입이 여러 그룹(CTRL, STIM)에 존재할 때, 
# 그 세포 타입에서 공통적으로 발현되는 마커 유전자를 찾고 싶을 때 사용
# nk.markers: STIM 그룹에서 NK 세포들의 공통 마커 유전자
nk.markers <- FindConservedMarkers(ifnb, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
head(nk.markers)
# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
# 각 세포 유형별로 대표 마커 유전자들이 어떻게 발현되는지를, CTRL vs STIM 조건 별로 나눠서 DotPlot 형태로 시각화
# factor(..., levels=...) : y축 (세포 유형 이름을 원하는 순서로 정렬)
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
    "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))
# markers.to.plot : x축 (시각화할 마커 유전자 목록 지정_각 유전자는 특정 세포 타입의 마커 유전자임)
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
# DotPlot으로 유전자 발현 시각화
# split.by = "stim" : 조건(ctrl vs stim)을 나눠서 그리라는 의미
DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()


# 조건에 따라 다르게 발현되는 유전자 식별
# Identify differential expressed genes across conditions
# 통합된 데이터를 기반으로, 자극 전후 세포 동일한 세포 유형의 유전자 발현 변화를 비교
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
# pseudobulk 만들기: 같은 세포 타입(e.g. "CD14 Mono") + 같은 조건(e.g. CTRL) 그룹에 속하는 모든 세포의 유전자 발현 평균
aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", "stim"), return.seurat = TRUE)
# 비교할 유전자 지정 : Interferon 자극에 반응하는 대표 유전자들
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
# CD14 Mono 자극 전후 발현량 산점도 => 자극(STIM) 후 강한 발현 증가 (ISG 계열 유전자 발현 증가)
p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
# CD4 Naive T 자극 전후 발현량 산점도 => 비슷한 유전자들이 반응하지만 변화 폭은 작음
p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
p2 + p4
# x=y: 변화X, x<y: STIM에서 발현량 높음, x>y: STIM에서 발현량 낮음
# => interferon 자극에 공통적으로 반응하지만, CD14 Mono가 더 강하게 반응.

# 새로운 메타데이터 컬럼인 celltype.stim 열 생성 
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
# "B_STIM"과 "B_CTRL" 집단 사이의 유전자 발현 차이(DEG) 계산
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)
FeaturePlot(ifnb, features = c("CD3D", "CXCL10", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey", "red"), reduction = "umap")
plots <- VlnPlot(ifnb, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "seurat_annotations", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


# SCTransform-normalized 데이터 세트로 통합 수행
# Perfom integration with SCTransform-normalized datasets
ifnb <- LoadData("ifnb")
# split datasets and process without integration
# ifnb[['RNA']] : ifnb 객체 안에 있는 RNA 레이어 즉, 원래의 유전자 발현값
# split(..., f = ifnb$stim) : stim조건(CTRL, STIM)에 따라 분할
# => RNA 발현 값을 자극 조건(STIM vs CTRL)에 따라 나눠서 저장
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunUMAP(ifnb, dims = 1:30)
# 비통합
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
# integrate datasets
# 하나의 RNA 레이어로 되어있던 것을 분할하여 CTRL/STIM으로 나눠 두 개의 RNA 레이어로 존재.
# IntegrateLayers() : 이 분리된 레이어를 보고, 서로 대응되는 세포 상태를 찾고, 공통 구조를 학습해서 조건 간 batch effect를 제거하면서 비교 가능한 축으로 표현
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, normalization.method = "SCT", verbose = F)
ifnb <- FindNeighbors(ifnb, reduction = "integrated.dr", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.6)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.dr")
# 통합
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
# perform differential expression
ifnb <- PrepSCTFindMarkers(ifnb)
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)