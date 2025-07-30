# Installation
library(devtools)
install.packages("remotes")  # 없으면 먼저 설치
install_github("navinlabcode/copykat")


# Prepare the readcount input file
library(Matrix)
library(data.table)
data_dir <- "C:/Users/user/Desktop/gbm_data/GBM2_extracted/GBM2_filtered_feature_bc_matrix"

# 1. matrix.mtx.gz 압축 해제 (임시 파일 위치 지정)
mtx_path <- file.path(data_dir, "matrix.mtx")
if (!file.exists(mtx_path)) {
  R.utils::gunzip(filename = file.path(data_dir, "matrix.mtx.gz"), destname = mtx_path, remove = FALSE)
}

# 2. sparse matrix 읽기
expr_matrix <- readMM(mtx_path)

# 3. gene 이름 읽기
features <- fread(cmd = paste("gzip -dc", file.path(data_dir, "features.tsv.gz")), header = FALSE)
gene_names <- features$V2  # 2열: gene symbol
# 4. cell barcode 읽기
barcodes <- fread(cmd = paste("gzip -dc", file.path(data_dir, "barcodes.tsv.gz")), header = FALSE)
cell_names <- barcodes$V1

# 5. 행렬에 rownames, colnames 지정
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- cell_names

# 6. dense matrix로 변환 (copyKAT 입력용)
expr_mat_dense <- as.matrix(expr_matrix)

# # 7. 확인
dim(expr_mat_dense)
head(rownames(expr_mat_dense))
head(colnames(expr_mat_dense))
# expr_mat_dense[1:5, 1:5]

# 8. Seurat 오브젝트 생성
library(Seurat)
# 중복 유전자 제거
expr_mat_dense <- expr_mat_dense[!duplicated(rownames(expr_mat_dense)), ]
raw <- CreateSeuratObject(counts = expr_mat_dense, project = "copykat.test", min.cells = 0, min.features = 0)

# 9. 원본 raw count matrix 얻기 (counts 슬롯)
exp.rawdata.gbm <- as.matrix(GetAssayData(raw, assay = "RNA", slot = "counts"))

# 10. 셀 2000개 무작위 샘플링
set.seed(42)
sampled.cells <- sample(colnames(exp.rawdata.gbm), 2000)
# 11. 샘플 셀로 subset 생성
exp.subset <- exp.rawdata.gbm[, sampled.cells]
# # 12. subset 행렬을 파일로 저장
# write.table(
#   exp.subset,
#   file = "C:/Users/user/Desktop/CopyKAT/exp.rawdata.subset.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = TRUE,
#   col.names = TRUE
# )


# Running copykat
library(copykat)
copykat.test <- copykat(
    rawmat=exp.subset, 
    id.type="S", 
    ngene.chr=3, 
    win.size=10, 
    KS.cut=0.1, 
    sam.name="test", 
    distance="euclidean", 
    norm.cell.names="",
    output.seg="FALSE", 
    plot.genes="TRUE", 
    genome="hg20",
    n.cores=1
)

pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")),]  ##keep defined cells
CNA.test <- data.frame(copykat.test$CNAmat)
# 결과 저장
library(data.table)
fwrite(pred.test, "C:/Users/user/Desktop/CopyKAT/gbm2_output/pred_test.csv")
fwrite(CNA.test, "C:/Users/user/Desktop/CopyKAT/gbm2_output/CNA_test.csv")


# navigate prediction results
head(pred.test)
#                            cell.names copykat.pred
# TGCGGGTCATTGTGCA.1 TGCGGGTCATTGTGCA-1      diploid
# ATCCGAACAAGGTTCT.1 ATCCGAACAAGGTTCT-1    aneuploid
# CGTTGGGGTTTAGCTG.1 CGTTGGGGTTTAGCTG-1    aneuploid
# TACAGTGTCAACACCA.1 TACAGTGTCAACACCA-1    aneuploid
# ACTATCTTCTCACATT.1 ACTATCTTCTCACATT-1      diploid
# GTGAAGGAGATCTGCT.1 GTGAAGGAGATCTGCT-1      diploid
head(CNA.test[ , 1:5])
#   chrom chrompos  abspos TGCGGGTCATTGTGCA.1 ATCCGAACAAGGTTCT.1
# 1     1  1042457 1042457        -0.05135258        -0.01826672
# 2     1  1265484 1265484        -0.05135258        -0.01826672
# 3     1  1519859 1519859        -0.05135258        -0.01826672
# 4     1  1826619 1826619        -0.05135258        -0.01826672
# 5     1  2058465 2058465        -0.05135258        -0.01826672
# 6     1  2280372 2280372        -0.05135258        -0.01826672


# heatmap 생성
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")


# (Optional) define subpopulations of aneuploid tumor cells
tumor.cells <- rownames(pred.test)[pred.test$copykat.pred == "aneuploid"]
tumor.mat <- CNA.test[, colnames(CNA.test) %in% tumor.cells]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')


# UMAP
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Seurat 객체 생성
seu <- CreateSeuratObject(counts = exp.subset)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:20)

# 2. 메타데이터 추가
# pred.test: copykat.pred 정보 담긴 데이터프레임
# > head(colnames(seu))
# [1] "115_2_70"    "124_1_341"   "143_7d_557"  "105_B1_2561" "105A_1176"  
# [6] "102_1457"
# > head(rownames(pred.test))
# [1] "X115_2_70"   "X124_1_341"  "X143_7d_557" "X105A_1176"  "X102_1457"  
# [6] "X125_2_311"
# 1. pred.test rownames의 접미사 ".1"을 "-1"로 변환
rownames(pred.test) <- sub("\\.1$", "-1", rownames(pred.test))
common_cells <- intersect(colnames(seu), rownames(pred.test))
seu <- subset(seu, cells = common_cells)
# 2. copykat.pred 메타데이터 추가
seu$copykat.pred <- pred.test[common_cells, "copykat.pred"]
# 3. tumor.cells 정의 (aneuploid 세포명)
tumor.cells <- rownames(pred.test)[pred.test$copykat.pred == "aneuploid"]
# 4. hc.umap(서브클론 클러스터링 결과) 가정, clones.pred 생성 및 이름 맞추기
clones.pred <- hc.umap
names(clones.pred) <- tumor.cells
# 5. seu에 subclone 컬럼 초기화 및 할당
seu$subclone <- NA
common_cells_subclone <- intersect(colnames(seu), names(clones.pred))
seu$subclone[match(common_cells_subclone, colnames(seu))] <- clones.pred[common_cells_subclone]
# 6. factor 변환 (시각화 편의)
seu$subclone <- factor(seu$subclone)
# 7. 시각화 - 전체 세포 copykat.pred로
DimPlot(seu, group.by = "copykat.pred", pt.size = 0.6, label = TRUE) + ggtitle("CNV status")
# 8. 시각화 - 종양 세포만 subclone으로
tumor_seu <- subset(seu, subset = copykat.pred == "aneuploid" & !is.na(subclone))
DimPlot(tumor_seu, group.by = "subclone", pt.size = 0.6, label = TRUE) + ggtitle("Tumor Subclones")


# aneuploid와 diploid 간의 DEG 분석
# 종양세포(aneuploid)와 정상세포(diploid)의 유전자 발현 차이 파악
# 어떤 유전자들이 종양세포에서 특이적으로 높거나 낮게 발현되는지 알 수 있음.
# 종양세포의 특징적인 분자 마커(molecular markers)를 찾을 수 있음.
library(Seurat)
Idents(seu) <- seu$copykat.pred
deg_results <- FindMarkers(seu,
                           ident.1 = "aneuploid",
                           ident.2 = "diploid",
                           min.pct = 0.25,
                           logfc.threshold = 0.25)
head(deg_results)
# > head(deg_results)
#                p_val avg_log2FC pct.1 pct.2     p_val_adj
# OSMR   9.810405e-134   2.843446 0.894 0.051 3.589725e-129
# C1R    2.631813e-133   3.836900 0.862 0.048 9.630066e-129
# FKBP10 1.144237e-127   3.738053 0.780 0.024 4.186877e-123
# NNMT   5.400547e-125   3.490604 0.915 0.087 1.976114e-120
# CALD1  2.383748e-122   3.198034 0.935 0.112 8.722372e-118
# MEG3   9.282073e-122   2.873868 0.833 0.052 3.396403e-117
# 결과를 csv로 저장 (원하면)
write.csv(deg_results, file = "DEG_aneuploid_vs_diploid.csv")
# 상위 5개 유전자 feature plot
top5_genes <- rownames(head(deg_results, 5))
FeaturePlot(seu, features = top5_genes, min.cutoff = "q10", max.cutoff = "q90", cols = c("lightgrey", "red"))