# Treatment condition에 따른 beta-cell 내부 잠재 substructure 탐색

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse

# 필요한 info가 포함된 dataset 만들기
adata = sc.read_h5ad('C:/Users/user/Desktop/T1D/GSE211799_adata_atlas.h5ad')

# cell type & condition 정의
# beta cell 필터
is_beta = adata.obs['cell_type_integrated_v2_parsed'] == 'beta'
# T2D treatment 조건
treatment_conditions = [
    "T2D_mSTZ-treated_GLP-1_estrogen",
    "T2D_mSTZ-treated_estrogen",
    "T2D_mSTZ-treated_GLP-1",
    "T2D_mSTZ-treated_GLP-1_estrogen+insulin"
]
# 'condition' 컬럼 생성
adata.obs['condition'] = adata.obs['CXG-DATA_diabetes_model']
# beta + treatment 조건 필터링
adata_beta_treat = adata[is_beta & adata.obs['condition'].isin(treatment_conditions)].copy()
#                           study_sample study file  reference  ...  CXG-DATA_sex_annotation cell_type_integrated_v2  cell_type_integrated_v2_parsed                                condition
# index                                                         ...
# CGGGTCATCATGGTCA-1-G8-STZ       STZ_G8   STZ   G8      False  ...             ground-truth                    beta                            beta  T2D_mSTZ-treated_GLP-1_estrogen+insulin
# CCTTACGCATGGAATA-1-G4-STZ       STZ_G4   STZ   G4      False  ...             ground-truth                    beta                            beta                   T2D_mSTZ-treated_GLP-1
# GTCATTTTCCGCGTTT-1-G6-STZ       STZ_G6   STZ   G6      False  ...             ground-truth                    beta                            beta          T2D_mSTZ-treated_GLP-1_estrogen
# GGAACTTGTGATAAGT-1-G4-STZ       STZ_G4   STZ   G4      False  ...             ground-truth                    beta                            beta                   T2D_mSTZ-treated_GLP-1
# GTAACTGCAATCTACG-1-G5-STZ       STZ_G5   STZ   G5      False  ...             ground-truth                    beta                            beta                T2D_mSTZ-treated_estrogen
# ...                                ...   ...  ...        ...  ...                      ...                     ...                             ...                                      ...
# GTTACAGTCGTCCGTT-1-G8-STZ       STZ_G8   STZ   G8      False  ...             ground-truth                    beta                            beta  T2D_mSTZ-treated_GLP-1_estrogen+insulin
# GGTGTTATCAGAGCTT-1-G4-STZ       STZ_G4   STZ   G4      False  ...             ground-truth                    beta                            beta                   T2D_mSTZ-treated_GLP-1
# CCAGCGACACACAGAG-1-G4-STZ       STZ_G4   STZ   G4      False  ...             ground-truth                    beta                            beta                   T2D_mSTZ-treated_GLP-1
# CGCTTCAAGATTACCC-1-G6-STZ       STZ_G6   STZ   G6      False  ...             ground-truth                    beta                            beta          T2D_mSTZ-treated_GLP-1_estrogen
# TAAACCGAGACGACGT-1-G8-STZ       STZ_G8   STZ   G8      False  ...             ground-truth                    beta                            beta  T2D_mSTZ-treated_GLP-1_estrogen+insulin
# [5834 rows x 56 columns]

# 확인
print(f"최종 beta cell 수: {adata_beta_treat.shape[0]}")   # 5834
print(adata_beta_treat.obs['condition'].value_counts())
# condition
# T2D_mSTZ-treated_GLP-1_estrogen            1670
# T2D_mSTZ-treated_estrogen                  1608
# T2D_mSTZ-treated_GLP-1_estrogen+insulin    1421
# T2D_mSTZ-treated_GLP-1                     1135
# Name: count, dtype: int64


# scICE preprocess
X = adata_beta_treat.X                  # gene expression matrix
gene_name = adata_beta_treat.var_names  # 유전자 이름 (index)
cell_name = adata_beta_treat.obs_names  # 세포 이름 (index)

# 2. 세포 필터링 기준 계산
min_genes_per_cell = 6
max_genes_per_cell = 0  # 0이면 제한 없음
mito_percent = 5.0
ribo_percent = 0.0
min_cells_per_gene = 15

# gene 이름이 Index이면 str으로 변환
gene_name = gene_name.to_series()

# --- 미토콘드리아 유전자 비율 계산 ---
bidx_mito = gene_name.str.upper().str.startswith("MT-").values
cell_mito_filter = (np.divide(X[:, bidx_mito].sum(axis=1).A1 if issparse(X) else np.sum(X[:, bidx_mito], axis=1), X.sum(axis=1).A1 if issparse(X) else np.sum(X, axis=1), out=np.zeros(X.shape[0]), where=(X.sum(axis=1).A1 if issparse(X) else np.sum(X, axis=1)) != 0) < mito_percent / 100) if mito_percent > 0 else np.ones(X.shape[0], dtype=bool)
# if mito_percent > 0:
#     mito_expr = X[:, bidx_mito].sum(axis=1).A1 if issparse(X) else np.sum(X[:, bidx_mito], axis=1)
#     total_expr = X.sum(axis=1).A1 if issparse(X) else np.sum(X, axis=1)
#     mito_ratio = np.divide(mito_expr, total_expr, out=np.zeros_like(mito_expr), where=total_expr != 0)
#     cell_mito_filter = mito_ratio < mito_percent / 100
# else:
#     cell_mito_filter = np.ones(X.shape[0], dtype=bool)

# --- 세포별 발현 유전자 수 계산 ---
gene_per_cell = np.diff(X.tocsr().indptr) if issparse(X) else np.count_nonzero(X, axis=1)
cell_gene_count_filter = ((gene_per_cell := (np.diff(X.tocsr().indptr) if issparse(X) else np.count_nonzero(X, axis=1))) >= min_genes_per_cell) & ((gene_per_cell < max_genes_per_cell) if max_genes_per_cell > 0 else True)
# cell_gene_count_filter = gene_per_cell >= min_genes_per_cell
# if max_genes_per_cell > 0:
#     cell_gene_count_filter &= gene_per_cell < max_genes_per_cell

# --- 최종 cell 필터 ---
cell_filter = cell_mito_filter & cell_gene_count_filter
X = X[cell_filter]
cell_name = cell_name[cell_filter]

# 3. 유전자 필터링 기준 계산
X_csc = X.tocsc() if issparse(X) else X
cell_per_gene = np.diff(X_csc.indptr) if issparse(X) else np.count_nonzero(X, axis=0)
gene_filter = cell_per_gene >= min_cells_per_gene
X = X[:, gene_filter]
gene_name = gene_name[gene_filter]

# 4. 발현량 0인 유전자 제거
gene_sum = X.sum(axis=0).A1 if issparse(X) else np.sum(X, axis=0)
nonzero_gene_filter = gene_sum != 0
X = X[:, nonzero_gene_filter]
gene_name = gene_name[nonzero_gene_filter]

# 5. 유전자 평균 발현량 기준 정렬
gene_means = X.mean(axis=0).A1 if issparse(X) else np.mean(X, axis=0)
sort_idx = np.argsort(gene_means)
X = X[:, sort_idx]
gene_name = gene_name.iloc[sort_idx]

# 6. 밀집 행렬로 변환 및 DataFrame 생성
X_dense = X.toarray() if issparse(X) else X
pre_df = pd.DataFrame(X_dense, columns=gene_name.values)
pre_df.insert(0, 'cell', cell_name.values)

# 결과 확인
print(f"scICE input pre_df shape: {pre_df.shape}")   # (5834, 14663)
print(f"sparsity: {1 - (pre_df.iloc[:,1:] != 0).sum().sum() / pre_df.iloc[:,1:].size:.4f}")   # 0.7749

# 저장
pre_df.to_csv('C:/Users/user/Desktop/T1D/beta/t2d_beta_treatment_preprocess.csv', index=False)
