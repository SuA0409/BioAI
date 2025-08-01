# T2D 치료 조건 간 Beta-cell 반응 프로그램 차이 추출
# pip install cnmf
# pip install harmonypy
# pip install scikit-misc

from cnmf import cNMF
import numpy as np
import scanpy as sc

# cNMF 실행

cnmf_obj = cNMF(output_dir="C:/Users/user/Desktop/T1D/t2d_treatment/cnmf_output", name="t2d_cnmf")

cnmf_obj.prepare(
    counts_fn="C:/Users/user/Desktop/T1D/t2d_treatment/cnmf_data.txt",
    components=np.arange(5, 12),  # K=5~13 범위 탐색 # ori: (5, 14)
    n_iter=100,
    seed=42
)

cnmf_obj.factorize(worker_i=0, total_workers=1)
# total_workers: 데이터를 나누어 병렬 연산할 때 전체 몇 개의 작업으로 나눌지
# worker_i: 현재 이 프로세스가 몇 번째 워커인지 (0부터 시작)
# -> 병렬로 여러 개를 돌리면 속도는 빨라지지만, RAM/CPU 리소스를 많이 먹음.

cnmf_obj.combine()

cnmf_obj.k_selection_plot()
# 최적 K 선정용 plot 출력

cnmf_obj.consensus(k=9, density_threshold=0.01)
# 최적의 K값으로 consensus clustering 수행

usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(K=9, density_threshold=0.01)