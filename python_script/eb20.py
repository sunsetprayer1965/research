# eb_heteroN.py  —— 每列用自己的 N_k（推荐）
import numpy as np

# 百分比 -> 比例
perc = np.array([
    8, 8.5, 9.333333333333334, 8.25, 9, 9, 9, 8.625, 8.444444444444445, 8.4,
    8.09090909090909, 8.166666666666666, 8.307692307692308, 8.428571428571429,
    8.4, 8.625, 8.588235294117647, 8.666666666666668, 9.05263157894737, 9.25
], float)
N = np.arange(100, 2001, 100, float)   # 各自的 N_k

r = perc / 100.0
Y = r * N

r_bar = r.mean()
S2 = r.var(ddof=1)
mean_invN = np.mean(1.0 / N)

tau2_hat = max(S2 - r_bar * (1.0 - r_bar) * mean_invN, 0.0)
if tau2_hat == 0.0:
    v_hat = 1e12  # 近似全收缩
else:
    v_hat = max((r_bar * (1.0 - r_bar)) / tau2_hat - 1.0, 1e-6)

alpha_hat = r_bar * v_hat
beta_hat  = (1.0 - r_bar) * v_hat

p_eb = (alpha_hat + Y) / (alpha_hat + beta_hat + N)
print("\t".join(f"{x*100:.6f}" for x in p_eb))
