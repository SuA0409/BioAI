import gzip
import pandas as pd
import random
import numpy as np

# 파일 경로
filepath = 'C:/Users/user/Desktop/Z4eq.csv.gz'

# Step 1: 헤더만 읽어서 열 이름 추출
with gzip.open(filepath, 'rt') as f:
    header_line = f.readline().strip()
    all_columns = header_line.split(',')

# Step 2: 'cell' 열 제외하고 랜덤으로 100개 선택
gene_columns = all_columns[1:]  # 첫 번째는 cell name이라고 가정
random_genes = random.sample(gene_columns, 100)
selected_columns = ['cell'] + random_genes  # cell 열 포함

# Step 3: usecols에 전달해서 원하는 열만 로딩
df = pd.read_csv(filepath, usecols=selected_columns)


## preprocess

# 기본 정보 추출
cell_name = df.cell.values  # cell type
# print(cell_name.head(5))  # 'cell' 열의 내용을 출력
gene_name = df.columns[1:]  # 유전자 이름
# print(gene_name)  # 2번째 열부터의 내용을 출력

# 데이터 크기/희소도 확인 
size_df = df.shape  # 데이터프레임의 크기 (3994, 101)
sp_ndf = 1.0 - (df.iloc[:, 1:] != 0).sum().sum() / df.iloc[:, 1:].size  # 희소도 계산 : 0.9659
print(f"data size: {size_df}, sparsity: {sp_ndf:.4f}")

# 희소 행렬/밀집 행렬 선택
from scipy.sparse import csc_matrix, issparse
X = df.iloc[:, 1:].to_numpy(dtype=np.float32) if sp_ndf < 0.3 else csc_matrix(df.iloc[:, 1:].values.astype(np.float32))
# sp_ndf가 0.3보다 작아서 csr_matrix로 변환 : csr_matrix(희소행렬)은 0이 많을수록 효율적으로 저장할 수 있는 구조
# print(X)
# print(f'각 행 기준 누적합: {X.indptr}')  # csr_matrix의 경우 각 행의 시작 인덱스
# print(f'열 위치 인덱스: {X.indices}')  # csr_matrix의 경우 각 비제로 요소의 열 인덱스
# print(f'저장된 데이터: {X.data}')  # csr_matrix의 경우 비제로 요소의 값

# 희소 행렬이면 df로 변환 (X: 희소행렬, df: 밀집행렬)
if issparse(X):
    df = pd.DataFrame(X.toarray(), columns=gene_name)
    # print(df.head(5))  # 변환된 데이터프레임의 첫 5행 출력
else:
    df = df.copy()  # 밀집 행렬이면 그대로 사용


# 필터링 조건 계산
# 함수 default 값
min_tp_c, min_tp_g=0, 0
max_tp_c, max_tp_g=float('Inf'), float('Inf')
min_genes_per_cell=6
# default값으로 200으로 설정되어있었지만, 현재 데이터로 적용했을 때 그 기준에 맞는 세포가 없었음
max_genes_per_cell=0
min_cells_per_gene=15
mito_percent=5.
ribo_percent=0

# 유전자별 nonzero cell 수 (세로합)
if issparse(X):
    # 희소행렬은 csc 형식으로 변환 후 계산 (열 axis=0 기준)
    X = X.tocsc()  # csc_matrix로 변환
    n_cell_counts = np.diff(X.indptr)  # 유전자별 nonzero cell(발현되는 세포) 수
    n_cell_counts_sum = X.sum(axis=0).A1  # 유전자별 총 발현량
    # A1은 1차원 배열로 변환
else:
    n_cell_counts = np.count_nonzero(X, axis=0)
    n_cell_counts_sum = np.sum(X, axis=0)
# 유전자 필터링 조건
bidx_1 = n_cell_counts_sum > min_tp_g # 유전자별 최소 발현량 (shape: (100,))
bidx_2 = n_cell_counts_sum < max_tp_g # 유전자별 최대 발현량 (shape: (100,))
bidx_3 = n_cell_counts >= min_cells_per_gene # 유전자별 최소 세포 수 (shape: (100,))
## 최종 gene index
fg_idx = bidx_1 & bidx_2 & bidx_3  # 최종 유전자 인덱스 # bidx_3와 동일한 출력값 (shape: (100,))

# 세포별 nonzero gene 수 (가로합)
if issparse(X):
    # 희소행렬은 csr 형식으로 변환 후 계산 (행 axis=1 기준)
    X = X.tocsr()
    n_gene_counts = np.diff(X.indptr)  # 세포별 nonzero gene(발현되는 유전자) 수
    n_gene_counts_sum = X.sum(axis=1).A1
else:
    n_gene_counts = np.count_nonzero(X, axis=1)
    n_gene_counts_sum = np.sum(X, axis=1)
# 세포별 필터링 조건
bidx_1 = n_gene_counts_sum > min_tp_c  # (shape: (3994,))
bidx_2 = n_gene_counts_sum < max_tp_c  # (shape: (3994,))
bidx_3 = n_gene_counts >= min_genes_per_cell  # (shape: (3994,))

# mitochochondrial gene filtering
bidx_mito = gene_name.str.contains(r'^(?i)mt-.')  
# "mt-"로 시작하고 그 뒤에 글자 최소 한개 있는 유전자 이름 (대소문자 구분X)
# mito_percent가 0이면 mito 유전자 필터링을 하지 않음 BUT, mito_percent = 5.0
if mito_percent == 0:
    bidx_4 = np.ones_like(bidx_1, dtype=bool)
else:
    mito_sum = X[:, bidx_mito].sum(axis=1).A1 if issparse(X) else np.sum(X[:, bidx_mito], axis=1)
    # 각 셀에서 미토콘드리아 관련 유전자 발현 합
    mito_total_sum = X.sum(axis=1).A1 if issparse(X) else np.sum(X, axis=1)
    # 각 셀에서 전체 유전자 발현 합
    # 나눗셈 방어 처리 (Julia에서는 연산에 오류가 발생하지 않아 따로 처리하지 않음)
    with np.errstate(divide='ignore', invalid='ignore'):
        mito_ratio = np.true_divide(mito_sum, mito_total_sum)
        # mito_otal_sum이 0인 경우가 아니라면 정상적으로 비율이 나옴.
        mito_ratio[~np.isfinite(mito_ratio)] = 0 
        # mito_total_sum이 0일 때만 NaN이나 inf가 생겨서 그 부분을 0으로 바꿔줌
    bidx_4 = mito_ratio < mito_percent / 100  # (shape: (3994,))

# ribosomal gene filtering
bidx_ribo = gene_name.str.contains(r'^(?i)RP[SL].')
# "RPS.." or "RPL.."로 시작하고 그 뒤에 글자가 최소 한개 있는 유전자 이름 (대소문자 구분X)
if ribo_percent == 0:
    bidx_5 = np.ones_like(bidx_1, dtype=bool)
else:
    ribo_sum = X[:, bidx_ribo].sum(axis=1).A1 if issparse(X) else np.sum(X[:, bidx_ribo], axis=1)
    ribo_total_sum = X.sum(axis=1).A1 if issparse(X) else np.sum(X, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        ribo_ratio = np.true_divide(ribo_sum, ribo_total_sum)
        ribo_ratio[~np.isfinite(ribo_ratio)] = 0 
    bidx_5 = ribo_ratio < ribo_percent / 100  # (shape: (3994,))

# max gene per cell filtering
if max_genes_per_cell == 0:
    bidx_6 = np.ones_like(bidx_1, dtype=bool)  # (shape: (3994,))
    # 최대 유전자 수 제한이 없으면, 모든 셀 통과(True)
else:
    bidx_6 = n_gene_counts < max_genes_per_cell
    # 최대 유전자 수 제한이 있으면, 각 셀의 발현 유전자 수가 기준보다 작은 셀만 True
## 최종 cell index
fc_idx = bidx_1 & bidx_2 & bidx_3 & bidx_4 & bidx_5 & bidx_6  # (shape: (3994,))


# 위의 QC 기준을 만족하는 cell과 gene만 필터링
oo_X = X[fc_idx][:, fg_idx]  # 필터링된 gene expression matrix (shape: (341, 64))
# min값을 200으로 설정했을 때, shape: (0, 64) : 데이터는 없음

# non-zero gene만 다시 필터링 즉, 하나라도 값이 없는 유전자 제거
nn_idx = (oo_X.sum(axis=0).A1 if issparse(oo_X) else np.sum(oo_X, axis=0)) != 0  # shape: (64,)
oo_X = oo_X[:, nn_idx]  # shape: (341, 64): 필터링 대상인 빈 열이 없어서 유지.
## 최종 필터링된 유전자 이름
norm_gene = gene_name[fg_idx][nn_idx]
# shape: (64,): fg_idx로 필터링된 100개 유전자 선택 후, 그 결과에 대해 nn_idx로 빈 열 제거하여 최종 유전자 이름 선택

# gene expression 평균 기준 정렬
gene_means = np.array(oo_X.mean(axis=0)).flatten() if issparse(oo_X) else np.mean(oo_X, axis=0)
# oo_X 값 형태 예시: (340, 47)  1.0
# 340: cell index(샘플 인덱스), 47: gene index(유전자 인덱스), 1.0: 그 cell에서 해당 유전자의 발현량
# 즉, 340번 cell에서 47번 유전자의 발현량이 1.0이라는 의미
# oo_X.mean(axis=0) : 각 유전자(ex. 47번째 유전자)마다 cell 개수(340개의 셀)에 대한 평균 발현량
# ex. 0번 유전자는 341개 셀에서 평균 발현값이 0.073, ...
s_idx = np.argsort(gene_means)  # argsort: gene_means 값을 오름차순 정렬했을 때의 인덱스 배열 반환
sorted_gene_names = norm_gene[s_idx]  # s_idx 순서대로 유전자 이름 정렬

## 최종 데이터프레임 정리
oo_X_dense = oo_X.toarray() if issparse(oo_X) else oo_X   # 희소 행렬 ->  밀집 행렬(array)
o_df = pd.DataFrame(oo_X_dense[:, s_idx], columns=sorted_gene_names)
# oo_X_dense[:, s_idx] : oo_X_dense에서 gene_means 기준으로 발현량 낮은 유전자 -> 높은 유전자 순으로 정렬
# columns=sorted_gene_names : 정렬된 유전자 이름을 열 이름으로 사용
# o_df: (cell x gene) 형태의 정렬된 발현값 테이블
o_df.insert(0, 'cell', cell_name[fc_idx])
# 첫 번째 열(index 0)에 'cell' 열 추가
# cell_name[fc_idx] : cell 이름 중, QC 필터링 통과한 셀 인덱스에 해당하는 cell 이름만 선택
# 출력확인
size_odf = o_df.shape  # shape: (341, 65)
sp_odf = 1.0 - (o_df.iloc[:, 1:] != 0).sum().sum() / o_df.iloc[:, 1:].size
# o_df에서 값이 0인 셀의 비율(희소도)

### main 함수 호출 : pre_df = scLENS.preprocess(ndf)
# preprocess(전처리) 후 반환되는 값은 o_df (pre_df = o_df)

# o_df.to_csv('C:/Users/user/Desktop/o_df.csv', index=False)