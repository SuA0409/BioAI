# T2D 치료 조건 간 Beta-cell 반응 프로그램 차이 추출
import numpy as np
import pandas as pd
import scanpy as sc

# 데이터 불러오기
adata = sc.read_h5ad('C:/Users/user/Desktop/T1D/GSE211799_adata_atlas.h5ad')
beta_cell_names = adata.obs_names[adata.obs['cell_type_integrated_v2_parsed'] == 'beta']
print(f"Found {len(beta_cell_names)} beta cells")
adata.obs['condition'] = adata.obs['CXG-DATA_diabetes_model']

# condition 컬럼 확인 및 분포 체크
print(f"Condition 분포: {adata.obs['condition'].value_counts()}")
# Condition 분포: condition
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

# T2D 저항성(대조군1), T2D(대조군2) T2D treatment(실험군들) condition만 추출
t2d_conditions = [
    "T2D_db/db", "T2D_mSTZ", "T2D_mSTZ-treated_estrogen", "T2D_mSTZ-treated_GLP-1_estrogen",
    "T2D_mSTZ-treated_GLP-1_estrogen+insulin", "T2D_mSTZ-treated_insulin", "T2D_mSTZ-treated_GLP-1"
    ]
adata_t2d = adata[adata.obs['condition'].isin(t2d_conditions)].copy()
print(adata_t2d.obs['condition'].value_counts())
# condition
# T2D_db/db                                  18230
# T2D_mSTZ-treated_GLP-1_estrogen             9863
# T2D_mSTZ                                    9184
# T2D_mSTZ-treated_estrogen                   7280
# T2D_mSTZ-treated_GLP-1                      6292
# T2D_mSTZ-treated_GLP-1_estrogen+insulin     4303
# T2D_mSTZ-treated_insulin                    4235
# Name: count, dtype: int64

# downsampling
# 조건1: 총 5000개 세포로 다운샘플링
# 조건2: Condition별 비율 유지
# 조건3: 각 condition 내에서 유전자 수 (n_genes_by_counts)가 
#       평균 ± 1.96 × 표준편차 내(신뢰구간 95%)인 세포만 샘플링 대상
# 조건4: 각 condition에서 해당 수만큼 무작위 추출

# 1. 조건 확인
total_cells = 5000  # 총 다운샘플링 대상 수

# 2. 조건별 cell 수와 비율 계산
condition_counts = adata_t2d.obs['condition'].value_counts()
condition_ratios = condition_counts / condition_counts.sum()
sample_sizes = (condition_ratios * total_cells).astype(int)

# 3. n_genes_by_counts 계산 (필요한 경우만)
if 'n_genes_by_counts' not in adata_t2d.obs.columns: adata_t2d.obs['n_genes_by_counts'] = (adata_t2d.X > 0).sum(axis=1).A1 if hasattr(adata_t2d.X, "A1") else (adata_t2d.X > 0).sum(axis=1)
# if 'n_genes_by_counts' not in adata_t2d.obs.columns:
#     if hasattr(adata_t2d.X, "A1"):
#         adata_t2d.obs['n_genes_by_counts'] = (adata_t2d.X > 0).sum(axis=1).A1
#     else:
#         adata_t2d.obs['n_genes_by_counts'] = (adata_t2d.X > 0).sum(axis=1)

# 4. 조건별 샘플링 (95% 신뢰구간 내에서)
filtered_indices = []
filtered_indices = sum([
    (in_ci.index.tolist() if len(in_ci := adata_t2d.obs[(adata_t2d.obs['condition'] == condition) & 
    ( (gene_counts := adata_t2d.obs[adata_t2d.obs['condition'] == condition]['n_genes_by_counts']).between(
        (mu := gene_counts.mean()) - 1.96 * (sigma := gene_counts.std()), 
        mu + 1.96 * sigma
    ))]) < n_sample 
    else np.random.choice(in_ci.index, size=n_sample, replace=False).tolist())
    for condition, n_sample in sample_sizes.items()
], [])
# for condition, n_sample in sample_sizes.items():
#     subset = adata_t2d.obs[adata_t2d.obs['condition'] == condition]
#     gene_counts = subset['n_genes_by_counts']

#     # 95% 신뢰구간 계산
#     mu, sigma = gene_counts.mean(), gene_counts.std()
#     lower, upper = mu - 1.96 * sigma, mu + 1.96 * sigma
#     in_ci = subset[(gene_counts >= lower) & (gene_counts <= upper)]

#     if len(in_ci) < n_sample:
#         print(f"{condition}: {n_sample}개 필요하지만, 신뢰구간 내 {len(in_ci)}개만 존재 → 전부 사용")
#         selected = in_ci.index
#     else:
#         selected = np.random.choice(in_ci.index, size=n_sample, replace=False)

#     filtered_indices.extend(selected)


# 5. 최종 다운샘플링된 AnnData 객체 생성
adata_t2d_sampled = adata_t2d[filtered_indices].copy()

# 6. 결과 확인
print(adata_t2d_sampled)
print(adata_t2d_sampled.obs['condition'].value_counts())

# 7. 저장
adata_t2d_sampled.write("C:/Users/user/Desktop/T1D/t2d_treatment_5000.h5ad")