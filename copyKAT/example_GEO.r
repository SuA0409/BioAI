# Installation
library(devtools)
install.packages("remotes")  # 없으면 먼저 설치
install_github("navinlabcode/copykat")
# install.packages("copykat")

# Prepare the readcount input file
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928
# 1. 압축 풀기
untar("C:/Users/user/Desktop/CopyKAT/GSE131928_RAW.tar", exdir = "C:/Users/user/Desktop/CopyKAT/unpacked/")
# 2. TSV 파일 불러오기
library(data.table)
tsv_path <- "C:/Users/user/Desktop/CopyKAT/unpacked/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz"
expr_mat <- data.table::fread(cmd = paste("gzip -dc", tsv_path), data.table = FALSE)
rownames(expr_mat) <- expr_mat[,1]
expr_mat <- expr_mat[,-1]
expr_mat <- as.matrix(expr_mat)
# 3. Seurat 오브젝트 생성
library(Seurat)
raw <- CreateSeuratObject(counts = expr_mat, project = "copykat.test", min.cells = 0, min.features = 0)
exp.rawdata <- as.matrix(GetAssayData(raw, assay = "RNA", slot = "counts"))
# 4. 파일 저장
set.seed(42)  # 재현 가능성
sampled.cells <- sample(colnames(exp.rawdata), 500)  # 셀 500개 무작위 추출
exp.subset <- exp.rawdata[, sampled.cells]  # subset 생성
# 지정된 형식으로 저장 (tab-separated, quote 없음, rownames 포함)
# write.table(
#   exp.subset,
#   file = "C:/Users/user/Desktop/CopyKAT/exp.rawdata.subset.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = TRUE
# )

# Running copykat
library(copykat)
copykat.test <- copykat(rawmat=exp.subset, id.type="S", ngene.chr=3, win.size=10, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FALSE", plot.genes="TRUE", genome="hg20",n.cores=1)

pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")),]  ##keep defined cells
CNA.test <- data.frame(copykat.test$CNAmat)
# local에 저장
library(data.table)
fwrite(pred.test, "pred_test.csv")
fwrite(CNA.test, "CNA_test.csv")


# navigate prediction results
head(pred.test)
#             cell.names copykat.pred
# X115_2_70     115_2_70      diploid
# X124_1_341   124_1_341    aneuploid
# X143_7d_557 143_7d_557    aneuploid
# X105A_1176   105A_1176      diploid
# X102_1457     102_1457    aneuploid
# X125_2_311   125_2_311    aneuploid
head(CNA.test[ , 1:5])
#   chrom chrompos  abspos   X115_2_70    X124_1_341
# 1     1  1042457 1042457  0.39214822  0.5223714002
# 2     1  1265484 1265484  0.39214822  0.5223714002
# 3     1  1519859 1519859  0.25832599  0.4213610275
# 4     1  1826619 1826619  0.12450375  0.3203506548
# 5     1  2058465 2058465  0.12450375  0.3203506548
# 6     1  2280372 2280372 -0.05056901 -0.0001013315

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
# 1. pred.test rownames 앞 'X' 제거 및 seu subset
rownames(pred.test) <- sub("^X", "", rownames(pred.test))
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
# 결과를 csv로 저장 (원하면)
write.csv(deg_results, file = "DEG_aneuploid_vs_diploid.csv")
# 상위 5개 유전자 feature plot
top5_genes <- rownames(head(deg_results, 5))
FeaturePlot(seu, features = top5_genes, min.cutoff = "q10", max.cutoff = "q90", cols = c("lightgrey", "red"))


