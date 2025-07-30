####### Data visualization methods in Seurat #######
install.packages("C:/Users/user/Desktop/Seurat/tutorial/3/pbmc3k.SeuratData_3.1.4.tar.gz", repos=NULL, type="source")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
# RNA assay를 포함하는 Seurat 객체 생성됨 (2700개의 PBMC 단일세포 포함)
pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")
# 2700개(세포 수)의 세포 각각에 "group1" 또는 "group2"를 랜덤으로 하나씩 배정. ("groups"라는 새 column 생성됨)
pbmc3k.final$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k.final), replace = TRUE)
# 분석에 사용할 특정 유전자(PBMC에서 주요하게 발현되는 마커 유전자들) 목록을 벡터로 저장
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final


# Five visualizations of marker feature expression
# Feature의 발현 패턴 시각화
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)
# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(pbmc3k.final, features = features)
# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(pbmc3k.final, features = features)
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(pbmc3k.final, features = features) + RotatedAxis()
# Single cell heatmap of feature expression
DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)


# New additions to FeaturePlot
# Plot a legend to map colors to expression levels
FeaturePlot(pbmc3k.final, features = "MS4A1")
# Adjust the contrast in the plot
FeaturePlot(pbmc3k.final, features = "MS4A1", min.cutoff = 1, max.cutoff = 3)
# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
FeaturePlot(pbmc3k.final, features = c("MS4A1", "PTPRCAP"), min.cutoff = "q10", max.cutoff = "q90")
# Visualize co-expression of two features simultaneously
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), blend = TRUE)
# Split visualization to view expression by groups (replaces FeatureHeatmap)
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), split.by = "groups")


# Updated and expanded visualization functions
# Violin plots can also be split on some variable. 
# Simply add the splitting variable to object metadata and pass it to the split.by argument
VlnPlot(pbmc3k.final, features = "percent.mt", split.by = "groups")
# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()
# DimPlot replaces TSNEPlot, PCAPlot, etc. 
# In addition, it will plot either 'umap', 'tsne', or 'pca' by default, in that order
DimPlot(pbmc3k.final)

pbmc3k.final.no.umap <- pbmc3k.final
pbmc3k.final.no.umap[["umap"]] <- NULL
DimPlot(pbmc3k.final.no.umap) + RotatedAxis()
# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. 
# This can be changed with the `group.by` parameter
DoHeatmap(pbmc3k.final, features = VariableFeatures(pbmc3k.final)[1:100], cells = 1:500, size = 4, angle = 90) + NoLegend()


# Applying themes to plots
baseplot <- DimPlot(pbmc3k.final, reduction = "umap")
# Add custom labels and titles
baseplot + labs(title = "Clustering of 2,700 PBMCs")
# Use community-created themes, overwriting the default Seurat-applied theme Install ggmin with remotes::install_github('sjessa/ggmin')
baseplot + ggmin::theme_powerpoint()   # ggmin 패키지 존재하지 않아서 설치 X -> 실행안됨
# Seurat also provides several built-in themes, such as DarkTheme; for more details see
# ?SeuratTheme
baseplot + DarkTheme()
# Chain themes together
baseplot + FontSize(x.title = 20, y.title = 20) + NoLegend()


# Interactive plotting features
# Include additional data to display alongside cell names by passing in a data frame of information.
# Works well when using FetchData
plot <- FeaturePlot(pbmc3k.final, features = "MS4A1")
# 마우스 커서를 셀(점) 위에 올리면, 해당 셀의 정보가 보여짐
# 이때, 셀마다 표시할 추가 정보를 지정할 수 있음 (ident: 각 셀의 cluster 식별자, PC_1: PCA의 첫 번째 주성분 값, nFeature_RNA: 각 셀의 유전자 feature 수)
HoverLocator(plot = plot, information = FetchData(pbmc3k.final, vars = c("ident", "PC_1", "nFeature_RNA")))

# cell cluster 중 "DC"라는 cluster 이름을 "CD14+ Mono"로 변경
pbmc3k.final <- RenameIdents(pbmc3k.final, DC = "CD14+ Mono")
plot <- DimPlot(pbmc3k.final, reduction = "umap")
# CellSelector()는 해당 plot 위에서 마우스로 영역을 클릭/드래그하여 셀을 선택
# select.cells에 선택된 셀들의 셀ID 벡터가 저장됨.
select.cells <- CellSelector(plot = plot)
head(select.cells)   # 마우스로 영역을 선택해야 select.cells에 셀ID가 저장됨.
Idents(pbmc3k.final, cells = select.cells) <- "NewCells"
# Now, we find markers that are specific to the new cells, and find clear DC markers
# FindMarkers()는 두 그룹 사이의 유전자 발현 차이 비교해서 특이적으로 발현된 유전자(marker)를 찾는 함수
# "NewCells"와 "CD14+ Mono" 그룹을 비교
newcells.markers <- FindMarkers(pbmc3k.final, ident.1 = "NewCells", ident.2 = "CD14+ Mono", min.diff.pct = 0.3, only.pos = TRUE)
# min.diff.pct = 0.3 : 두 그룹 사이 발현 비율 차이가 30% 이상인 유전자만 고려
# only.pos = TRUE : "NewCells"에서 더 많이 발현된 유전자만 출력
# => 즉, NewCells 그룹에서 특이적으로 높게 발현된 유전자들만 뽑음.
head(newcells.markers)


# Plotting Accessories
# LabelClusters and LabelPoints will label clusters (a coloring variable) or individual points on a ggplot2-based scatter plot
plot <- DimPlot(pbmc3k.final, reduction = "pca") + NoLegend()
# LabelClusters() : cluster별로 라벨링
LabelClusters(plot = plot, id = "ident")
# Both functions support `repel`, which will intelligently stagger labels and draw connecting lines from the labels to the points or clusters
# LabelPoints() : 셀들의 ID(e.g. "AAAATTGA-1")로 라벨링, repel = TRUE : 라벨이 겹치지 않도록 자동으로 간격을 두고 선으로 연결선 그려줌
LabelPoints(plot = plot, points = TopCells(object = pbmc3k.final[["pca"]]), repel = TRUE)

plot1 <- DimPlot(pbmc3k.final)
# Create scatter plot with the Pearson correlation value as the title
plot2 <- FeatureScatter(pbmc3k.final, feature1 = "LYZ", feature2 = "CCL5")
# Combine two plots
plot1 + plot2
# Remove the legend from all plots
(plot1 + plot2) & NoLegend()
