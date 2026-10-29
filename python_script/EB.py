#   u_hat = mean(Y^k / N)
#   v_hat = (N-1)/N * [u_hat*(1-u_hat) / var_hat(Y^k/N)] - 1
#   alpha_hat = u_hat * v_hat
#   beta_hat  = (1 - u_hat) * v_hat
#   p_eb^k = (alpha_hat + Y^k) / (alpha_hat + beta_hat + N)

import numpy as np
from pathlib import Path

# ====== Config ======
FILE_MDS  = "output_chroms_msplMDS_linking_number.txt"
FILE_NMDS = "output_chroms_msplNMDS_linking_number.txt"
FILE_PM1  = "output_chroms_msplPM1_linking_number.txt"
FILE_PM2  = "output_chroms_msplPM2_linking_number.txt"

N = 1000  # given

OUT_FILE = "eb_estimates.txt"

# ====== Helper ======
def load_file(path: str) -> np.ndarray:
    arr = np.loadtxt(path, dtype=float)
    arr = np.atleast_1d(arr).astype(float).ravel()
    arr = arr / 100.0

    if len(arr) != 120:
        raise ValueError(f"{path} 应有120个数字，实际 {len(arr)}")
    return arr

def eb_for_pair(r: np.ndarray, N: int) -> np.ndarray:
    """r: shape (4,), return EB (4,)"""
    u_hat = np.mean(r)
    print(f"u_hat: {u_hat}")
    var_r = np.var(r, ddof=1) if len(r) > 1 else 0.0
    print(f"var_r: {var_r}")
    v_hat = (N - 1) / N * (u_hat * (1 - u_hat) / var_r) - 1

    a_hat = u_hat * v_hat
    b_hat = (1 - u_hat) * v_hat
    Y = r * N
    p_eb = (a_hat + Y) / (a_hat + b_hat + N)
    return p_eb

# ====== Main ======
r_mds  = load_file(FILE_MDS)
r_nmds = load_file(FILE_NMDS)
r_pm1  = load_file(FILE_PM1)
r_pm2  = load_file(FILE_PM2)

eb_matrix = np.zeros((120, 4))
for i in range(120):
    r = np.array([r_mds[i], r_nmds[i], r_pm1[i], r_pm2[i]], dtype=float)
    eb_matrix[i, :] = eb_for_pair(r, N)
eb_matrix = eb_matrix*100.0


# 保存成4列txt
np.savetxt(OUT_FILE, eb_matrix, fmt="%.6f", delimiter="\t",
           header="MDS\tNMDS\tPM1\tPM2", comments='')

print(f"[完成] 已输出 EB 估计：{OUT_FILE}")

