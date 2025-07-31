# 1. T1D vs T2D 간 Beta-cell 유전자 발현 차이
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gget
import seaborn as sns

# 1. 데이터 불러오기
adata = sc.read("C:/Users/user/Desktop/T1D/beta_cells_condition_analysis.h5ad")
print(adata)

# 2. condition 컬럼 확인 및 분포 체크
adata = adata[adata.obs['condition'].isin(['T1D_NOD', 'T2D_db/db'])].copy()
print(f"Condition 분포: {adata.obs['condition'].value_counts()}")

# 3. 정규화 및 로그 변환
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw = adata   # 원본 보존

# 4. PCA & UMAP
# 원래는 주성분별 분산 비율 확인 후, n_comps 설정.
sc.pp.pca(adata, n_comps=20)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)   # 이웃 그래프 재계산
sc.tl.umap(adata, min_dist=0.2, spread=1.0)   # UMAP 재계산
sc.pl.umap(adata, color='condition', palette='Set1', title='UMAP (tuned): T1D vs T2D (Beta cells)')
# hyperparameter 
# n_comps : 6, 10, 20 
# n_neighbors: 10, 15, 20, 30
# n_pcs: 6, 8, 10, 15, 20
# min_dist: 0.2, 0.3, 0.5

# # 주성분별 분산 비율 확인 (explained variance ratio)
# explained_variance_ratio = adata.uns['pca']['variance_ratio']
# # 분산 비율 시각화
# plt.plot(range(1, len(explained_variance_ratio)+1), explained_variance_ratio, 'o-')
# plt.xlabel('Principal Component')
# plt.ylabel('Explained Variance Ratio')
# plt.title('PCA Explained Variance Ratio')
# plt.show()

# 5. clustering (beta-cell 내 이질성 탐색)
sc.tl.leiden(adata, resolution=0.5)
# 데이터 내의 세포들을 유사성에 따라 그룹으로 나눔(resolution 값이 클수록 더 많은 클러스터로 분리됨) 
sc.pl.umap(adata, color=['condition', 'leiden'], title=['Condition', 'Leiden Clusters'])
# >>> adata.obs['leiden']
# index
# CTAACTTGTGAGTGAC-1-SRR7610303-NOD_elimination    2
# TGGATGTTCTCGTCAC-1-MUC13639-VSG                  0
# CAGTAACAGGACATTA-1-SRR7610298-NOD_elimination    3
# CATAAGCGTGCCCAGT-1-MUC13639-VSG                  0
# GGCAATTCAGACGCCT-1-SRR7610298-NOD_elimination    3
#                                                 ..
# AGACTCACAGTGCCTG-1-MUC13641-VSG                  1
# GGGCTACAGGTTGAGC-1-MUC13639-VSG                  1
# CAGATCAAGTCGATAA-1-SRR7610302-NOD_elimination    2
# TACAGTGCAGCCTTGG-1-SRR7610303-NOD_elimination    2
# GGGACCTAGAAATTGC-1-MUC13641-VSG                  0
# Name: leiden, Length: 422, dtype: category
# Categories (5, object): ['0', '1', '2', '3', '4']

# 6. 클러스터별 마커 유전자 추출
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groupby='leiden', n_genes=10, sharey=False)   # 시각화
marker_df = sc.get.rank_genes_groups_df(adata, group=None)
# sc.tl.rank_genes_groups()로 저장된 마커 유전자 분석 결과를 판다스 DataFrame 형태로 가져오는 함수
marker_df['gene_symbol'] = marker_df['names'].map(adata.var['gene_symbol'])
print(marker_df.head())
#   group               names     scores  logfoldchanges         pvals     pvals_adj gene_symbol
# 0     0  ENSMUSG00000024516  11.070742        0.468618  1.739552e-28  1.838474e-24      Sec11c
# 1     0  ENSMUSG00000005873  10.566929        0.637660  4.241591e-26  3.362097e-22       Reep5
# 2     0  ENSMUSG00000003380  10.283041        0.480701  8.403355e-25  4.440613e-21      Rabac1
# 3     0  ENSMUSG00000092486  10.071440        1.740572  7.388860e-24  2.928390e-20         NaN
# 4     0  ENSMUSG00000031919  10.029382        1.637411  1.132229e-23  3.605808e-20       Tmed6

# 클러스터별 Top 5 마커 유전자 출력
[print(f"Cluster {g}:\n", marker_df.query("group == @g").nlargest(5, "logfoldchanges")[['gene_symbol', 'logfoldchanges', 'pvals_adj']], "\n") for g in marker_df['group'].unique()]


# 7. DEG 분석: T1D_NOD vs T2D_db/db within beta cells
sc.tl.rank_genes_groups(
    adata,
    groupby='condition',
    method='wilcoxon',
    groups=['T1D_NOD'],
    reference='T2D_db/db',
    use_raw=True,
    pts=True
)

# 8. DEG 시각화 및 저장
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, title='T1D vs T2D (Beta cells)')
deg_df = sc.get.rank_genes_groups_df(adata, group='T1D_NOD')
deg_df.to_csv("C:/Users/user/Desktop/T1D/t1d_vs_t2d_deg.csv", index=False)

# 9. 클러스터 라벨 시각화 (중심 좌표에 메인 cell type annotate)
umap_coords = adata.obsm['X_umap']
clusters = adata.obs['leiden'].cat.categories
fig, ax = plt.subplots(figsize=(6, 6))
sc.pl.umap(adata, color='leiden', ax=ax, show=False)

[ax.text(np.mean(adata.obsm['X_umap'][adata.obs['leiden'] == c, 0]),
         np.mean(adata.obsm['X_umap'][adata.obs['leiden'] == c, 1]),
         f'Cluster {c}', fontsize=10, weight='bold',
         ha='center', va='center',
         bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')) for c in adata.obs['leiden'].cat.categories]

ax.set_title('Beta Cell State Clusters')
plt.show()


# T1D, T2D 각각에서 클러스터 분포 비교하기
# 1.1) 조건별(leiden 클러스터 분포)
ct1 = pd.crosstab(adata.obs['condition'], adata.obs['leiden'], normalize='index')
print(ct1)
# leiden            0         1        2        3         4
# condition
# T1D_NOD    0.000000  0.000000  0.54321  0.45679  0.000000
# T2D_db/db  0.476923  0.457692  0.00000  0.00000  0.065385
# 1.2) 시각화 (비율 막대 그래프)
ct1.plot(kind='bar', stacked=True, figsize=(8,5), colormap='tab20')
plt.ylabel('Proportion')
plt.title('Leiden Cluster Distribution by Condition (T1D vs T2D)')
plt.legend(title='Leiden Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

# 2.1) 클러스터별 condition 비율 계산
ct2 = pd.crosstab(adata.obs['leiden'], adata.obs['condition'], normalize='index')
print(ct2)
# condition  T1D_NOD  T2D_db/db
# leiden
# 0              0.0        1.0
# 1              0.0        1.0
# 2              1.0        0.0
# 3              1.0        0.0
# 4              0.0        1.0
# 2.2) 막대그래프로 시각화 (스택 바)
ct2.plot(kind='bar', stacked=True, figsize=(8,5), colormap='Set1')
plt.ylabel('Proportion')
plt.xlabel('Leiden Cluster')
plt.title('Condition Distribution per Cluster (T1D vs T2D)')
plt.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()


# 클러스터별 marker gene 기반 상태 라벨링하기
# 1) 클러스터별 top 마커 유전자 확인
marker_df = sc.get.rank_genes_groups_df(adata, group=None)
marker_df['gene_symbol'] = marker_df['names'].map(adata.var['gene_symbol'])
print(marker_df.head())
# 클러스터별 top 5 유전자 출력
from IPython.display import display
# for cluster in marker_df['group'].unique():
#     print(f"Cluster {cluster} top 5 markers:")
#     display(marker_df.query('group == @cluster').nlargest(5, 'logfoldchanges')[['gene_symbol', 'logfoldchanges', 'pvals_adj']])
#     print('\n')
result_str = "\n".join([
    f"Cluster {cluster} top 5 markers:\n" +
    marker_df.query('group == @cluster').nlargest(5, 'logfoldchanges')[['gene_symbol', 'logfoldchanges', 'pvals_adj']].to_string(index=False) + "\n"
    for cluster in marker_df['group'].unique()
])
print(result_str)
with open("C:/Users/user/Desktop/T1D/cluster_top5_markers.txt", "w") as f:
    f.write(result_str)

# 2) 직접 상태 이름 붙이기 (수동!! 직접 해야함!!!)
# 클러스터 ID → 상태 이름 수동 매핑 (예시)
cluster_to_state = {
    '0': 'Mature Beta',
    '1': 'ER Stress Beta',
    '2': 'Proliferative Beta',
    '3': 'Dysfunctional Beta',
    '4': 'Immature Beta'
}
adata.obs['beta_state'] = adata.obs['leiden'].map(cluster_to_state).astype('category')

# 3) 라벨 시각화
sc.pl.umap(adata, color=['beta_state', 'condition'], title=['Beta Cell States', 'Condition'])