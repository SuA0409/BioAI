# package installation & library loading
install.packages("Matrix")
install.packages("Seurat")
install.packages("patchwork")
install.packages("BiocManager")
BiocManager::install("harmony")

library("Seurat")
library("dplyr")
library("patchwork")
library("harmony")
library("ggplot2")


# Data loading
HC_filt.matrix <- Read10X("C:/Users/user/Desktop/Seurat/0714_BK강의/HC1/")
HC <- CreateSeuratObject(counts = HC_filt.matrix, project = "HCA")

P25_filt.matrix <- Read10X("C:/Users/user/Desktop/Seurat/0714_BK강의/P25_T6/")
P25 <- CreateSeuratObject(counts = P25_filt.matrix, project = "HCA")


# 기본 QC
HC[["percent.mt"]] <- PercentageFeatureSet(HC, pattern = "^MT-")
VlnPlot(HC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot 결과 보고 아래 적합한 조건 설정
HC <- subset(HC, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 4)
                          # live cell          # live cell           # mt 높을수록 세포가 죽은 것

P25[["percent.mt"]] <- PercentageFeatureSet(P25, pattern = "^MT-")
VlnPlot(P25, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot 결과 보고 아래 적합한 조건 설정
P25<- subset(P25, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & percent.mt < 5)


# Dimension Reduction
HC <- NormalizeData(HC)
HC <- FindVariableFeatures(HC, selection.method = "vst", nfeatures = 4000, verbose = FALSE)
HC <- ScaleData(HC)
HC <- RunPCA(HC)
# Elbow plot to determine the number of PCs to use
ElbowPlot(HC)   # 분산이 급격히 감소하는 지점
# 3부터 이지만, 안정성을 위해 PC=5

P25 <- NormalizeData(P25)
P25 <- FindVariableFeatures(P25, selection.method = "vst", nfeatures = 4000, verbose = FALSE)
P25 <- ScaleData(P25)
P25 <- RunPCA(P25)
ElbowPlot(P25)


# UMAP 그리기
HC <- RunUMAP(HC, reduction = "pca", dims = 1:5)
HC <- FindNeighbors(HC, reduction = "pca", dims = 1:5)
HC <- FindClusters(HC)

P25 <- RunUMAP(P25, reduction = "pca", dims = 1:9)
P25 <- FindNeighbors(P25, reduction = "pca", dims = 1:9)
P25 <- FindClusters(P25)
# 웬만하면 dims를 10이상으로 잡기

DimPlot(HC, reduction = "umap",  label = TRUE)
DimPlot(P25, reduction = "umap", label = TRUE)


# Meta.data에 정보 추가하기기
HC$GEM <- "GEM1"
HC$Disease <- "HC"
str(HC@meta.data$Disease)

P25$GEM <- "GEM2"
P25$Disease <- "Sepsis"


# Integration 전처리

HC[["RNA"]] <- as(object = HC[["RNA"]], Class = "Assay5")
P25[["RNA"]] <- as(object = P25[["RNA"]], Class = "Assay5")

HC$orig.ident <- "GEM1"
P25$orig.ident <- "GEM2"

PBMC_combined <- merge(HC,
  y = list(P25),
  add.cell.id = c("HC", "P25"))
str(PBMC_combined)
str(PBMC_combined@meta.data)
str(PBMC_combined@assays$nCount_RNA)

# scale.data 제거
PBMC_combined@assays$RNA@layers$scale.data <- NULL
str(PBMC_combined)
str(PBMC_combined@assays$RNA@layers)
PBMC_combined <- JoinLayers(PBMC_combined)

# GEM 기준으로 RNA assay 분할
PBMC_combined[["RNA"]] <- split(PBMC_combined[["RNA"]], f = PBMC_combined$GEM)  # 중요
str(PBMC_combined)
# Normalize, HVG, Scale, PCA
PBMC_combined <- NormalizeData(PBMC_combined)
PBMC_combined <- FindVariableFeatures(PBMC_combined)
PBMC_combined <- ScaleData(PBMC_combined)
PBMC_combined <- RunPCA(PBMC_combined)
str(PBMC_combined)


# Integration, 함수에 F1 눌러보기
PBMC_combined <- IntegrateLayers(
  object = PBMC_combined, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)   # 중요


# Integrated UMAP 그리기
PBMC_combined <- JoinLayers(PBMC_combined)   # 중요


PBMC_combined <- FindNeighbors(PBMC_combined, reduction = "harmony", dims = 1:30)
PBMC_combined <- FindClusters(PBMC_combined, resolution = 0.8, cluster.name = "harmony_clusters")
PBMC_combined <- RunUMAP(PBMC_combined, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")

PBMC_combined <- RenameIdents(PBMC_combined, "0" = "Monocyte")
PBMC_combined <- RenameIdents(PBMC_combined, "1" = "T/NK")
PBMC_combined <- RenameIdents(PBMC_combined, "2" = "T/NK")
PBMC_combined <- RenameIdents(PBMC_combined, "3" = "Platelet")
PBMC_combined <- RenameIdents(PBMC_combined, "4" = "Platelet")
PBMC_combined <- RenameIdents(PBMC_combined, "5" = "RBC")
PBMC_combined <- RenameIdents(PBMC_combined, "6" = "T/NK")
PBMC_combined <- RenameIdents(PBMC_combined, "7" = "B cell")
PBMC_combined <- RenameIdents(PBMC_combined, "8" = "Platelet")
PBMC_combined <- RenameIdents(PBMC_combined, "9" = "RBC")
PBMC_combined <- RenameIdents(PBMC_combined, "10" = "RBC")
PBMC_combined <- RenameIdents(PBMC_combined, "11" = "Monocyte")
PBMC_combined <- RenameIdents(PBMC_combined, "12" = "T/NK")
PBMC_combined <- RenameIdents(PBMC_combined, "13" = "B cell")

DimPlot(PBMC_combined, reduction = "umap_harmony", label = TRUE)

View(PBMC_combined@meta.data)
# Annotation을 위한 Marker gene set 선정
str(PBMC_combined@active.ident)
genes_to_plot <- c("FAM138A", "AL669831.5", "OR4F16")
str(genes_to_plot)
FeaturePlot(PBMC_combined, features = c(
    "IGHM", "CD79A", "MS4A1",         # B
    "NKG7", "CD3D", "CD8A", "CD4",    # T/NK
    "CTSS", "LYZ", "SERPINA1",        # Monocyte
    "HBD", "HBM", "AHSP",            # RBC
    "PPBP", "PF4", "NRGN"         # Platelet
    ))

DotPlot(
  PBMC_combined,
  features = c(
  "IGHM", "CD79A", "MS4A1",         # B
  "NKG7", "CD3D", "CD8A", "CD4",    # T/NK
  "CTSS", "LYZ", "SERPINA1",        # Monocyte
  "PPBP", "PF4", "NRGN",            # RBC
  "HBD", "HBM", "AHSP"         # Platelet
  ))


# Annotation을 위한 Dot plot 찍어보기
DotPlot(
  PBMC_combined,
  features = c(genes_to_plot),
  group.by = "harmony_clusters"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10)
  ) +
  coord_flip()


# Annotation 확인하기 - Feature plot
p <- FeaturePlot(PBMC_combined, features = c(
  "IGHM", "CD79A", "MS4A1",         # B
  "NKG7", "CD3D", "CD8A", "CD4",    # T/NK
  "CTSS", "LYZ", "SERPINA1",        # Monocyte
  "PPBP", "PF4", "NRGN",            # RBC
  "HBD", "HBM", "AHSP"         # Platelet
  ))
str(p)

# PDF로 저장
# ggsave(filename = 빈칸1, plot = 빈칸2, width = 12, height = 10)


# Annotation 확인하기 - Vln plot
VlnPlot(
  PBMC_combined,
  features = c("CD79A", "CD3D", "SERPINA1", "PPBP", "HBD"),
  pt.size = 0.5,              # 점 크기 (0이면 점 생략)
  ncol = 3                  # 한 줄에 몇 개의 플롯을 배치할지
)

table(T_NK_Sepsis$Disease)
unique(T_NK_Sepsis$Disease)
# DEG 확인하기
str(PBMC_combined@active.ident)
T_NK_Sepsis <- subset(PBMC_combined, idents = "T/NK")  # Monocyte만 추출
str(T_NK_Sepsis)
View(T_NK_Sepsis@meta.data)   # 옆 팝업 창으로 T_NK_Sepsis의 meta.data 확인
# Active identity 설정
Idents(T_NK_Sepsis) <- "Disease"  # Disease를 기준으로 identity 설정
View(T_NK_Sepsis@active.ident)
# Sepsis vs HC 비교 (Sepsis가 case, HC가 control)
##F1 눌러서 함수 확인하기
T_NK_DEG <- FindMarkers(
  object = T_NK_Sepsis,
  ident.1 = "Sepsis",  # case 그룹
  ident.2 = "HC",      # control 그룹
  logfc.threshold = 0,  # 모든 DEG 확인 (필요 시 조절)
  min.pct = 0.1,        # 최소 발현 비율 (필요 시 조절)
  return.thresh = 0.  # 유의한 DEG만 필터링 (기본값은 0.05)
)

# 결과 확인
head(T_NK_DEG)


# DEG 결과를 CSV로 저장
write.csv(
  오브젝트,
  file = "빈칸",
  row.names = TRUE)