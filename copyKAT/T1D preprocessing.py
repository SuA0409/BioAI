import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse
# import gzip
# with gzip.open("C:/Users/user/Desktop/GSE211799_adata_atlas.h5ad.gz", 'rb') as f:
#     adata = ad.read_h5ad(f)

# Raw data 구조 및 구성 파악
adata = sc.read_h5ad('C:/Users/user/Desktop/T1D/GSE211799_adata_atlas.h5ad')

print(adata.shape)                      # (n_cells, n_genes) = (301796, 31706)
print(adata)                            # adata 주요 구성 내용
print(adata.obs.head())                 # 각 셀의 메타데이터: sample_id, condition, diabetes_status, etc.
print(adata.var.head())                 # 유전자 정보
print(adata.uns.keys())                 # 시각화용 설정값, 필드 설명 등
print(adata.obsm.keys())                # PCA, UMAP, latent space 등 low-dimensional embedding
print(adata.X[:5, :5])                  # 처음 5x5 행렬 출력
print(adata.uns["field_descriptions"])  # 메타데이터 항목 설명


# 원하는 정보만 저장
# 1. CopyKAT 용의 raw count(adata.X) + gene(adata.var_names) + cell info(adata.obs_names)
# raw count matrix 추출
raw_X = adata.raw.to_adata().X if adata.raw else adata.X
raw_genes = adata.raw.var_names if adata.raw else adata.var_names
raw_cells = adata.obs_names
# 가공된 AnnData 객체 생성
adata_copykat = ad.AnnData(X=raw_X)
adata_copykat.var_names = raw_genes
adata_copykat.obs_names = raw_cells
# 저장
np.random.seed(42)
sampled_cells = np.random.choice(adata_copykat.obs_names, size=5000, replace=False)
adata_copykat_subset = adata_copykat[sampled_cells].copy()
adata_copykat_subset.write("C:/Users/user/Desktop/T1D/copykat_input.h5ad")
# adata_copykat.write("copykat_input.h5ad")

# 2. DEG 분석용의 log-normalized expression + label 정보
# log-normalized expression matrix (이미 계산되어 있다고 가정)
deg_X = adata.X
deg_genes = adata.var_names.copy()
deg_cells = adata.obs_names.copy()
# 군집 또는 cell type 라벨
deg_labels = adata.obs['cell_type_integrated_v2_parsed'].copy()
# 가공된 AnnData 객체 생성
adata_deg = ad.AnnData(X=deg_X, obs=pd.DataFrame({'cell_type': deg_labels}, index=deg_cells))
adata_deg.var_names = deg_genes
adata_deg.obs_names = deg_cells
# 저장
np.random.seed(42)
sampled_cells = np.random.choice(adata_deg.obs_names, size=5000, replace=False)
adata_deg_subset = adata_deg[sampled_cells].copy()
adata_deg_subset.write("C:/Users/user/Desktop/T1D/deg_input.h5ad")
# adata_deg.write("deg_input.h5ad")

# 3. UMAP용 expression + UMAP 좌표
# expression + UMAP embedding
# 1. UMAP 좌표가 있는지 확인
umap_coord = adata.obsm['X_umap'] if 'X_umap' in adata.obsm else None
if umap_coord is None:
    raise ValueError("UMAP coordinates ('X_umap') not found in adata.obsm")
# 2. 셀 이름 5000개 샘플링
np.random.seed(42)
sampled_cells = np.random.choice(adata.obs_names, size=5000, replace=False)
# 3. subset만 추출
adata_subset = adata[sampled_cells].copy()
# 4. X를 float32로 변환 (희소 or dense 모두 대응)
adata_subset.X = adata_subset.X.astype("float32") if issparse(adata_subset.X) else np.asarray(adata_subset.X, dtype="float32")
# 5. UMAP 좌표 subset만 추출
umap_subset_coord = umap_coord[np.isin(adata.obs_names, sampled_cells), :]
# 6. 새 AnnData 생성
adata_vis_subset = ad.AnnData(
    X=adata_subset.X,
    obs=adata_subset.obs.copy(),
    var=adata.var.copy(),
    obsm={"X_umap": umap_subset_coord}
)
adata_vis_subset.var_names = adata.var_names
adata_vis_subset.obs_names = sampled_cells
# 7. 저장
adata_vis_subset.write("C:/Users/user/Desktop/T1D/umap_input.h5ad")
# adata_vis.write("umap_input.h5ad")

# 4. beta-cell subset
print(adata.obs['cell_type_integrated_v2_parsed'].value_counts())
# cell_type_integrated_v2_parsed
# beta             102143
# alpha             40935
# immune            31703
# E non-endo.       29177
# delta             24775
# stellate a.       18332
# endothelial       13469
# ductal             8742
# E endo.            7748
# gamma              6999
# beta+delta         5104
# stellate q.        4970
# alpha+delta        1901
# beta+gamma         1209
# delta+gamma        1069
# endo. prolif.       887
# lowQ                853
# alpha+beta          683
# schwann             617
# acinar              480
# Name: count, dtype: int64
# 1. beta cell 이름만 먼저 추출
beta_cell_names = adata.obs_names[adata.obs['cell_type_integrated_v2_parsed'] == 'beta']
print(f"Found {len(beta_cell_names)} beta cells")
# 2. 샘플링 크기 결정
n_sample = 5000 if len(beta_cell_names) >= 5000 else len(beta_cell_names)
# 3. 무작위 샘플링
np.random.seed(42)
sampled_cells = np.random.choice(beta_cell_names, size=n_sample, replace=False)

# 4. 샘플 셀로 subset & 복사 (메모리 절감)
beta_cells_subset = adata[sampled_cells].copy()
beta_cells_subset.write("C:/Users/user/Desktop/T1D/beta_cells_subtype_analysis.h5ad")

# 5. beta-cell subset + T1D/control 조건 라벨 정보
# 조건 컬럼 이름 예시: 'CXG-DATA_diabetes_model' (적절히 바꿔 써야 함)
# 5. condition 컬럼 추가
beta_cells_subset.obs['condition'] = beta_cells_subset.obs['CXG-DATA_diabetes_model']
beta_cells_subset.write("C:/Users/user/Desktop/T1D/beta_cells_condition_analysis.h5ad")


# 확인
# 1. CopyKAT 용의 raw count(adata.X) + gene(adata.var_names) + cell info(adata.obs_names)
adata_copykat = sc.read("C:/Users/user/Desktop/T1D/copykat_input.h5ad")
print(adata_copykat)   # n_obs × n_vars = 5000 × 31706
print(adata_copykat.obs_names[:5])      # 셀 이름 확인
# Index(['GGGAGATAGAGCAATT-1-SRR7610301-NOD_elimination',
#        'AGTGGGAGTAACGTTC-1-E14_5-embryo',
#        'AAAGCAAGTGCTGTAT-1-SRR10751515-spikein_drug',
#        'AAATGCCGTGTTGGGA-1-G8-STZ', 'CACCAGGAGTTGTCGT-1-E13_5-embryo'],
#       dtype='object', name='index')
print(adata_copykat.var_names[:5])      # 유전자 이름 확인
# Index(['ENSMUSG00000000001', 'ENSMUSG00000000003', 'ENSMUSG00000000028',
#        'ENSMUSG00000000031', 'ENSMUSG00000000037'],
#       dtype='object', name='EID')
print(type(adata_copykat.X), adata_copykat.X.shape)  # raw count matrix
# <class 'scipy.sparse._csr.csr_matrix'> (5000, 31706)

# 2. DEG 분석용의 log-normalized expression + label 정보
adata_deg = sc.read("C:/Users/user/Desktop/T1D/deg_input.h5ad")
print(adata_deg)   # n_obs × n_vars = 5000 × 31706, obs: 'cell_type'
print(adata_deg.obs.head())             # 메타데이터 (cluster, label 등)
#                                                  cell_type
# index
# GGGAGATAGAGCAATT-1-SRR7610301-NOD_elimination        alpha
# AGTGGGAGTAACGTTC-1-E14_5-embryo                E non-endo.
# AAAGCAAGTGCTGTAT-1-SRR10751515-spikein_drug           beta
# AAATGCCGTGTTGGGA-1-G8-STZ                             beta
# CACCAGGAGTTGTCGT-1-E13_5-embryo                E non-endo.
print(adata_deg.var.head())             # 유전자 정보
# Empty DataFrame
# Columns: []
# Index: [ENSMUSG00000000001, ENSMUSG00000000003, ENSMUSG00000000028, ENSMUSG00000000031, ENSMUSG00000000037]
print(type(adata_deg.X), adata_deg.X.shape)  # log-normalized expression matrix
# <class 'scipy.sparse._csr.csr_matrix'> (5000, 31706)

# 3. UMAP용 expression + UMAP 좌표
adata_umap = sc.read("C:/Users/user/Desktop/T1D/umap_input.h5ad")
print(adata_umap)
print("UMAP 좌표 포함?", "X_umap" in adata_umap.obsm)
print(adata_umap.obsm["X_umap"][:5] if "X_umap" in adata_umap.obsm else "❌ 없음")
# [[ 6.1372433  -0.47747895]
#  [-2.8647733  -2.1728954 ]
#  [-3.3513782   2.934936  ]
#  [-5.7986574   4.881866  ]
#  [-4.7680893   4.7482347 ]]

# 4. beta-cell subset 
adata_beta = sc.read("C:/Users/user/Desktop/T1D/beta_cells_subtype_analysis.h5ad")
print(adata_beta)
print(adata_beta.obs.head())  # beta cell의 메타정보
#                                                              study_sample  ... cell_type_integrated_v2_parsed
# index                                                                      ...

# ACTTCCGCACACGCCA-1-MUC13640-VSG                              VSG_MUC13640  ...
#   beta
# ACGATGTAGATCCTAC-1-MUC13634-VSG                              VSG_MUC13634  ...
#   beta
# CTAACTTGTGAGTGAC-1-SRR7610303-NOD_elimination  NOD_elimination_SRR7610303  ...
#   beta
# CAACGATAGTCACTCA-1-MUC13640-VSG                              VSG_MUC13640  ...
#   beta
# TGGCGCAGTCCGAAGA-1-G1-STZ                                          STZ_G1  ...
#   beta
print(adata_beta.var.head())
#                    gene_symbol  used_integration  ... CXG-DATA_present_STZ  gene_symbol_FINAL
# EID                                               ...
# ENSMUSG00000000001       Gnai3             False  ...                 True              Gnai3
# ENSMUSG00000000003        Pbsn             False  ...                 True               Pbsn
# ENSMUSG00000000028       Cdc45             False  ...                 True              Cdc45
# ENSMUSG00000000031         H19             False  ...                 True                H19
# ENSMUSG00000000037       Scml2             False  ...                 True              Scml2
# 5. beta-cell subset + T1D/control 조건 라벨 정보
adata_beta_cond = sc.read("C:/Users/user/Desktop/T1D/beta_cells_condition_analysis.h5ad")
print(adata_beta_cond)
print("컬럼 목록:", adata_beta_cond.obs.columns.tolist())
print(adata_beta_cond.obs["condition"].value_counts())  # T1D vs control 분포
# condition
# T2D_db/db-treated_VSG                      474
# T2D_db/db-treated_PairFeed                 462
# T2D_db/db                                  260
# T1D_NOD                                    162
# T1D_NOD_prediabetic                         98
# T2D_mSTZ-treated_estrogen                   84
# T2D_mSTZ-treated_GLP-1_estrogen             83
# T2D_mSTZ-treated_GLP-1_estrogen+insulin     68
# T2D_mSTZ                                    66
# T2D_mSTZ-treated_insulin                    56
# T2D_mSTZ-treated_GLP-1                      55
# Name: count, dtype: int64

# expression 값을 DataFrame처럼 보기
df = pd.DataFrame(adata.X[:5, :5].todense() if issparse(adata.X) else adata.X[:5, :5],
                  index=adata.obs_names[:5], columns=adata.var_names[:5])
print(df)
# EID                                            ENSMUSG00000000001  ...  ENSMUSG00000000037
# index                                                              ...
# GGGAGATAGAGCAATT-1-SRR7610301-NOD_elimination                 0.0  ...                 0.0
# AGTGGGAGTAACGTTC-1-E14_5-embryo                               1.0  ...                 0.0
# AAAGCAAGTGCTGTAT-1-SRR10751515-spikein_drug                   3.0  ...                 0.0
# AAATGCCGTGTTGGGA-1-G8-STZ                                     0.0  ...                 0.0
# CACCAGGAGTTGTCGT-1-E13_5-embryo                               1.0  ...                 0.0