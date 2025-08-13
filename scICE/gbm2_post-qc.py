import os
import gzip
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

path = "C:/Users/user/Desktop/gbm_data/GBM2_extracted/GBM2_filtered_feature_bc_matrix"

# cell barcode 리스트 읽기
with gzip.open(os.path.join(path, "barcodes.tsv.gz"), "rt") as f:
    barcodes = [line.strip() for line in f]  # 11857

# gene 이름 리스트 읽기 (features.tsv.gz에서 2번째 컬럼)
with gzip.open(os.path.join(path, "features.tsv.gz"), "rt") as f:
    genes = [line.strip().split("\t")[1] for line in f]  # 36601

# matrix.mtx.gz 희소 행렬 읽기
mtx = mmread(os.path.join(path, "matrix.mtx.gz")).tocsr()  # CSR 희소 행렬로 변환  # (36601, 11857)

# 희소 행렬의 행, 열 차원 확인 (mtx.shape == (genes, cells) 인 경우 transpose 필요)
# 보통 10x mtx는 (genes x cells) 형태이므로, (cells x genes)로 전치 필요
print("Original matrix shape:", mtx.shape)  # (genes, cells)
mtx = mtx.transpose().tocsr()
print("Transposed matrix shape:", mtx.shape)  # (cells, genes)

# pandas DataFrame 생성 (큰 데이터는 메모리 문제 있으니 주의)
df_expr = pd.DataFrame.sparse.from_spmatrix(mtx, index=barcodes, columns=genes)
df_gbm = df_expr.reset_index().rename(columns={'index': 'cell'})
print(df_gbm.shape)  # (cell 수, gene 수)
print(df_gbm)

import pandas as pd
from zipfile import ZipFile
from io import StringIO, BytesIO

# 인덱스를 'cell'이라는 컬럼으로 옮기기
df_gbm = df_expr.reset_index().rename(columns={'index': 'cell'})

# 메모리상에 CSV 문자열 저장
csv_buffer = BytesIO()
df_gbm.to_csv(csv_buffer, index=False)

# ZIP 파일 생성 (CSV는 메모리에서 직접 삽입)
with ZipFile('gbm_data.zip', 'w') as zip_file:
    zip_file.writestr('gbm_data_pre.csv', csv_buffer.getvalue())

print("저장 완료: gbm_data.zip")


# import numpy as np
# import random

# # 랜덤 인덱스 5개 추출 (cell, gene 각각)
# cell_indices = random.sample(range(len(barcodes)), 5)
# gene_indices = random.sample(range(len(genes)), 5)

# # numpy 배열로 변환 (희소행렬은 toarray 또는 tocsr로 변환)
# mtx_dense = mtx.toarray()

# # mtx에서 선택된 위치 값 (5x5 행렬)
# vals_matrix = mtx_dense[np.ix_(cell_indices, gene_indices)]

# # df_expr에서도 같은 위치 값 추출
# cells = [barcodes[i] for i in cell_indices]
# genes_sel = [genes[j] for j in gene_indices]
# vals_df = df_expr.loc[cells, genes_sel].values

# # 두 배열 비교 (전체)
# print("Are values equal?:", np.array_equal(vals_matrix, vals_df))

# # 비교 결과 출력 (원하는 경우)
# print("Matrix values:\n", vals_matrix)
# print("DataFrame values:\n", vals_df)
