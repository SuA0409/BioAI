# T2D 치료 조건 간 Beta-cell 반응 프로그램 차이 추출
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1. 데이터 불러오기
adata = sc.read_h5ad('C:/Users/user/Desktop/T1D/t2d_treatment_5000.h5ad')

# 2. condition 컬럼 확인 및 분포 체크
print(f"Condition 분포: {adata.obs['condition'].value_counts()}")
# Condition 분포: condition
# T2D_db/db                                  1534
# T2D_mSTZ-treated_GLP-1_estrogen             830
# T2D_mSTZ                                    773
# T2D_mSTZ-treated_estrogen                   612
# T2D_mSTZ-treated_GLP-1                      529
# T2D_mSTZ-treated_GLP-1_estrogen+insulin     362
# T2D_mSTZ-treated_insulin                    356
# Name: count, dtype: int64

# 3. before filtering, analysis
# 세포별 발현 유전자 수 확인
# adata.obs에 세포별 발현 유전자 수 컬럼 추가 (0보다 큰 유전자 개수)
adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(axis=1).A1  # A1: sparse matrix → ndarray flatten
# 분포 시각화
plt.hist(adata.obs['n_genes_by_counts'], bins=50)
plt.xlabel('Number of genes expressed per cell')
plt.ylabel('Number of cells')
plt.title('Distribution of genes expressed per cell')
plt.show()   # 700 ~4700
# 유전자별 발현 세포 수 확인
# 각 유전자가 발현된 세포 수 (열 방향으로 0보다 큰 값 개수)
gene_cell_counts = (adata.X > 0).sum(axis=0).A1
plt.hist(gene_cell_counts, bins=50)
plt.xlabel('Number of cells expressing gene')
plt.ylabel('Number of genes')
plt.title('Distribution of cells expressing each gene')
plt.show()
# before filtering
print(f"전체 세포 수: {adata.n_obs}")   # 4996
print(f"전체 유전자 수: {adata.n_vars}")   # 31706

# 4. filtering
# 세포별 발현 유전자 수 계산 (adata 기준)
adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, "A1") else (adata.X > 0).sum(axis=1)
# 세포 필터링 (700 이상 4700 이하)
cells_before = adata.n_obs
adata = adata[(adata.obs['n_genes_by_counts'] >= 700) &
                (adata.obs['n_genes_by_counts'] <= 4700)].copy()
cells_after = adata.n_obs
print(f"세포 필터링: {cells_before} -> {cells_after} (700~4700 발현 유전자 수)")   # 4532
# 유전자별 발현 세포 수 계산
min_cells = max(int(cells_after * 0.01), 5)
genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=min_cells)
genes_after = adata.n_vars
print(f"유전자 필터링: {genes_before} -> {genes_after} (≥ {min_cells} cells)")   # 13460

# 5. after filtering, identifing 'n_genes_by_counts' 분포
print(adata.obs['n_genes_by_counts'].describe())
# print(adata.obs['condition'].value_counts())
# condition
# T2D_db/db                                  1239
# T2D_mSTZ-treated_GLP-1_estrogen             794
# T2D_mSTZ                                    746
# T2D_mSTZ-treated_estrogen                   588
# T2D_mSTZ-treated_GLP-1                      505
# T2D_mSTZ-treated_GLP-1_estrogen+insulin     343
# T2D_mSTZ-treated_insulin                    317
# Name: count, dtype: int64

# 6. 필터링된 데이터 저장
output_path = "C:/Users/user/Desktop/T1D/t2d_treatment_filtered.h5ad"
adata.write(output_path)

# adata = sc.read_h5ad('C:/Users/user/Desktop/T1D/t2d_treatment_filtered.h5ad')
# 7. Raw counts matrix 추출
counts_matrix = adata.raw.X if adata.raw is not None else adata.X   # gene x cell

# 8. counts matrix, cell, gene 이름 저장
# counts matrix를 text 파일로 저장 (cNMF에서 불러올 수 있도록)
# CSR sparse matrix → dense array
if hasattr(counts_matrix, "toarray"):
    counts_matrix = counts_matrix.toarray()
np.savetxt("C:/Users/user/Desktop/T1D/t2d_treatment/counts_matrix.txt", counts_matrix.T, fmt='%d', delimiter='\t')   # cell x gene
# gene names 저장
with open("C:/Users/user/Desktop/T1D/t2d_treatment/genes.txt", 'w') as f:
    for gene in adata.var_names:
        f.write(gene + "\n")
# cell names 저장
with open("C:/Users/user/Desktop/T1D/t2d_treatment/cells.txt", 'w') as f:
    for cell in adata.obs_names:
        f.write(cell + "\n")