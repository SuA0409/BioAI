import scanpy as sc
import zipfile
import os

# zip 파일 경로와 해제할 디렉토리 지정
zip_path = "C:/Users/user/Desktop/gbm_data/GBM2_filtered_feature_bc_matrix.zip"
extract_path = "C:/Users/user/Desktop/gbm_data/GBM2_extracted"

# zip 파일 해제
with zipfile.ZipFile(zip_path, 'r') as zip_ref:
    zip_ref.extractall(extract_path)

path = "C:/Users/user/Desktop/gbm_data/GBM2_extracted/GBM2_filtered_feature_bc_matrix/"
print(os.listdir(path))  # ['barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']

# Scanpy로 10X 형식 데이터 읽기
adata = sc.read_10x_mtx("C:/Users/user/Desktop/gbm_data/GBM2_extracted/GBM2_filtered_feature_bc_matrix/", var_names='gene_symbols')

# 결과 확인
print(adata)   # n_obs × n_vars = 11857 × 36601   var: 'gene_ids', 'feature_types'
print("Genes:", adata.var_names[:10].tolist())
print("Cell:", adata.obs_names[:10].tolist())


import gzip
import os

# 파일이 있는 폴더 경로
path = "C:/Users/user/Desktop/gbm_data/GBM2_extracted/GBM2_filtered_feature_bc_matrix/"

print("=== barcodes.tsv.gz ===")
with gzip.open(os.path.join(path, "barcodes.tsv.gz"), "rt") as f:
    for i in range(5):
        print(f"Cell {i+1}:", f.readline().strip())
# Cell 1: AAACCTGAGATGAGAG-1
# Cell 2: AAACCTGAGCCATCGC-1
# Cell 3: AAACCTGAGCTAGTGG-1
# Cell 4: AAACCTGAGGAATTAC-1
# Cell 5: AAACCTGAGGCAGGTT-1

print("\n=== features.tsv.gz ===")
with gzip.open(os.path.join(path, "features.tsv.gz"), "rt") as f:
    for i in range(5):
        print(f"Gene {i+1}:", f.readline().strip().split("\t"))
# Gene 1: ['ENSG00000243485', 'MIR1302-2HG', 'Gene Expression']
# Gene 2: ['ENSG00000237613', 'FAM138A', 'Gene Expression']
# Gene 3: ['ENSG00000186092', 'OR4F5', 'Gene Expression']
# Gene 4: ['ENSG00000238009', 'AL627309.1', 'Gene Expression']
# Gene 5: ['ENSG00000239945', 'AL627309.3', 'Gene Expression']

print("\n=== matrix.mtx.gz ===")
with gzip.open(os.path.join(path, "matrix.mtx.gz"), "rt") as f:
    for i in range(10):
        print(f"Line {i+1}:", f.readline().strip())
# Line 1: %%MatrixMarket matrix coordinate integer general
# Line 2: %metadata_json: {"software_version": "cellranger-7.1.0", "format_version": 2}
# Line 3: 36601 11857 28649568
# Line 4: 17 1 1
# Line 5: 25 1 2
# Line 6: 45 1 6
# Line 7: 46 1 3
# Line 8: 54 1 1
# Line 9: 61 1 2
# Line 10: 63 1 2
