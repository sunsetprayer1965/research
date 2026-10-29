import numpy as np
from itertools import combinations

# =========================================================
# 工具函数：边界检查与轻微裁剪（论文需 0<p<1，否则导数会发散）
# =========================================================
def _check_open01(x, name="p"):
    if np.any((x <= 0) | (x >= 1)):
        raise ValueError(
            f"{name} 含有不在开区间(0,1)内的值。论文假设 0<p<1。"
            f" 发现 min={x.min():.4g}, max={x.max():.4g}"
        )

def maybe_clip_open01(x, strict=True, eps=1e-4):
    """
    strict=True: 严格要求 0<p<1，遇到边界直接报错（论文原始假设）
    strict=False: 为工程稳健性，极轻微裁剪到 (eps,1-eps)
    推荐：strict=False 且 eps 至少不小于 1/(2N)
    """
    x = np.asarray(x, float)
    if strict:
        _check_open01(x)
        return x
    # 非严格：轻微裁剪
    low = eps
    high = 1.0 - eps
    return np.clip(x, low, high)

# =========================================================
# 变换及导数（arcsin sqrt 变换），论文中用于方差稳定化
# =========================================================
def arcsin_sqrt(p):
    # 仅在 0<p<1 定义良好
    return np.arcsin(np.sqrt(p))

def d_arcsin_sqrt_dp(p):
    # d/dp arcsin(sqrt(p)) = 1 / [2 sqrt(p(1-p))]
    return 1.0 / (2.0 * np.sqrt(p * (1.0 - p)))

# =========================================================
# 两类 G 函数（半软阈值化），以及它们对 u 的导数
# u = |z|/c
# =========================================================
def G1(u):
    # 连续平滑：小 u 收缩强，大 u 接近 1
    return (u*u) / (1.0 + u*u)

def dG1_du(u):
    return 2.0 * u / (1.0 + u*u)**2

def G2(u):
    # 带“硬”阈值：u<=0.5 则 0，u>0.5 时是平滑函数
    out = np.zeros_like(u)
    m = u > 0.5
    v = (u[m] - 0.5)
    out[m] = (v*v) / (1.0 + v*v)
    return out

def dG2_du(u):
    # 在 u<=0.5 处导数为 0（次梯度取 0 可行）
    out = np.zeros_like(u)
    m = u > 0.5
    v = (u[m] - 0.5)
    out[m] = 2.0 * v / (1.0 + v*v)**2
    return out

G_FUNCS = {"G1": (G1, dG1_du), "G2": (G2, dG2_du)}

# =========================================================
# 统计量：未变换 与 变换 两种路径
# 返回：
#   delta: 差值（未变换：p-q；变换：arcsin√p - arcsin√q）
#   z: 标准化差值（依据对应的方差近似）
#   ddelta_dp, ddelta_dq: delta 对 p/q 的导数
#   dz_dp, dz_dq: z 对 p/q 的导数（用于 Stein 校正）
#   及需要的 Var_p, Var_q（未变换时）
# 说明：此处按你要求假设 Np=Nq=N=常数=1000
# =========================================================
def stats_untransformed(p, q, N):
    """
    未变换路径：
        delta = p - q
        Var(delta) = p(1-p)/N + q(1-q)/N
        z = delta / sqrt(Var(delta))
    """
    delta = p - q
    Var = p*(1.0 - p)/N + q*(1.0 - q)/N
    S = np.sqrt(Var)
    z = delta / S

    # Var 对 p/q 的导数
    dVar_dp = (1.0 - 2.0*p) / N
    dVar_dq = (1.0 - 2.0*q) / N

    # z 对 p/q 的导数：z = delta / S
    dz_dp = (1.0 / S) - delta * (0.5 * dVar_dp) / (S**3)
    dz_dq = (-1.0 / S) - delta * (0.5 * dVar_dq) / (S**3)

    # delta 对 p/q 的导数
    ddelta_dp = np.ones_like(p)
    ddelta_dq = -np.ones_like(q)

    # 供 Stein 校正项使用
    Var_p = p*(1.0 - p)/N
    Var_q = q*(1.0 - q)/N

    return delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, Var_p, Var_q

def stats_transformed(p, q, N):
    """
    变换(arcsin√p)路径：
        tp = arcsin√p, tq = arcsin√q
        delta = tp - tq
        Var(tp) ≈ 1/(4N), Var(tq) ≈ 1/(4N)
        Var(delta) ≈ 1/(2N)
        z = delta / sqrt(1/(2N)) = delta * sqrt(2N)
    """
    tp = arcsin_sqrt(p)
    tq = arcsin_sqrt(q)
    delta = tp - tq
    z = delta * np.sqrt(2.0 * N)

    # d(delta)/dp, d(delta)/dq
    ddelta_dp = d_arcsin_sqrt_dp(p)
    ddelta_dq = -d_arcsin_sqrt_dp(q)

    # z = delta * sqrt(2N) ⇒ ∂z = sqrt(2N) * ∂delta
    dz_dp = np.sqrt(2.0 * N) * ddelta_dp
    dz_dq = np.sqrt(2.0 * N) * ddelta_dq

    return delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq

# =========================================================
# 收缩后的差值 δ̃ 及其导数
#   δ̃ = δ * G(u),  u = |z|/c
#   ∂δ̃/∂p = (∂δ/∂p)G(u) + δ * G'(u) * (1/c) * sgn(z) * ∂z/∂p
#   ∂δ̃/∂q 同理
# =========================================================
def delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du):
    u = np.abs(z) / c
    G_u = G(u)
    dG = dG_du(u)

    sign_z = np.sign(z)  # z=0 时取 0，连续延拓合理
    d_absz_dp = sign_z * dz_dp
    d_absz_dq = sign_z * dz_dq

    delta_tilde = delta * G_u

    ddt_dp = ddelta_dp * G_u + delta * dG * (1.0/c) * d_absz_dp
    ddt_dq = ddelta_dq * G_u + delta * dG * (1.0/c) * d_absz_dq
    return delta_tilde, ddt_dp, ddt_dq

# =========================================================
# 论文的 F̂(c) （两条路径）
# =========================================================
def Fhat_untransformed(p, q, N, c, G, dG_du, strict=True, eps=1e-4):
    # 读入比例，并根据 strict/eps 处理 0/1
    p = maybe_clip_open01(np.asarray(p, float), strict=strict, eps=eps)
    q = maybe_clip_open01(np.asarray(q, float), strict=strict, eps=eps)

    delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, Var_p, Var_q = stats_untransformed(p, q, N)
    delta_tilde, ddt_dp, ddt_dq = delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du)

    term1 = np.sum(delta_tilde**2)
    term2 = -2.0 * np.sum(delta * delta_tilde)
    # Stein 修正（未变换）
    term3 = 2.0 * np.sum(Var_p * ddt_dp - Var_q * ddt_dq)

    return float(term1 + term2 + term3), delta_tilde

def Fhat_transformed(p, q, N, c, G, dG_du, strict=True, eps=1e-4):
    p = maybe_clip_open01(np.asarray(p, float), strict=strict, eps=eps)
    q = maybe_clip_open01(np.asarray(q, float), strict=strict, eps=eps)

    delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq = stats_transformed(p, q, N)
    delta_tilde, ddt_dp, ddt_dq = delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du)

    term1 = np.sum(delta_tilde**2)
    term2 = -2.0 * np.sum(delta * delta_tilde)
    # 变换情形的 Stein 修正： (2/(4N)) * Σ(∂δ̃/∂p̂ - ∂δ̃/∂q̂) = (1/(2N)) * Σ(...)
    term3 = (1.0 / (2.0 * N)) * np.sum(ddt_dp - ddt_dq)

    return float(term1 + term2 + term3), delta_tilde

# =========================================================
# 在网格上最小化 F̂(c) 得到 ĉ；并按论文阈值 0.005 计算 zero_ratio
# =========================================================
def minimize_over_c(p, q, N, mode, G_name,
                    c_min=0.001, c_max=60.0, n_grid=1000,
                    strict=False, eps=1e-4, zero_const=0.005):
    """
    参数：
      mode: "untransformed" 或 "transformed"
      G_name: "G1" 或 "G2"
      zero_const: 论文建议的阈值常数（默认 0.005）
      strict=False: 默认用轻微裁剪，避免 0/1 导致导数发散
      eps: 非严格模式下的裁剪幅度（建议 >= 1/(2N)）
    """
    G, dG = G_FUNCS[G_name]
    cs = np.linspace(c_min, c_max, n_grid)

    best_c, best_F, best_delta_t = None, np.inf, None
    for c in cs:
        if mode == "untransformed":
            F, dt = Fhat_untransformed(p, q, N, c, G, dG, strict=strict, eps=eps)
        else:
            F, dt = Fhat_transformed(p, q, N, c, G, dG, strict=strict, eps=eps)
        if F < best_F:
            best_c, best_F, best_delta_t = c, F, dt

    # 论文中的“零项”判据：|δ̃| < const（const=0.005）
    zero_ratio = float(np.mean(np.abs(best_delta_t) < zero_const))
    return float(best_c), float(best_F), zero_ratio

# =========================================================
# 主流程：对给定的一组 reconstruction 两两比较，并输出四种组合结果
# =========================================================
def run_all_variants(R, names, N,
                     c_range=(0.01, 10.0), n_grid=1000,
                     strict=False, eps=None, zero_const=0.005):
    """
    R: list[np.ndarray]，每个是一个 reconstruction 的比例向量（0~1）
    names: list[str]，每个向量的名字
    N: 样本量（你的设定：所有条目 N_i ≡ N = 1000）
    eps: 非严格裁剪幅度，若 None 则取 max(1/(2N), 1e-4)
    """
    if eps is None:
        eps = max(1.0/(2.0*N), 1e-4)

    results = {}
    for (a, b) in combinations(range(len(R)), 2):
        p = R[a]; q = R[b]
        key = f"{names[a]} vs {names[b]}"
        results[key] = {}
        for mode in ["untransformed", "transformed"]:
            for G_name in ["G1", "G2"]:
                c_hat, Fmin, zr = minimize_over_c(
                    p, q, N, mode, G_name,
                    c_min=c_range[0], c_max=c_range[1], n_grid=n_grid,
                    strict=strict, eps=eps, zero_const=zero_const
                )
                results[key][f"{mode}-{G_name}"] = {
                    "c_hat": c_hat,
                    "Fhat_min": Fmin,
                    "zero_ratio": zr
                }
    return results

def pretty_print(results, zero_const=0.005):
    print("\n=== Semi-soft thresholding（论文一致，含 Stein 修正） ===")
    print(f"(零项判据：|δ̃| < {zero_const})")
    for pair, d in results.items():
        print(f"\nComparison: {pair}")
        for k, v in d.items():
            print(f"  {k:<20s}  c*={v['c_hat']:.3f}   Fhat_min={v['Fhat_min']:.6f}   zero_ratio={v['zero_ratio']:.2%}")

# =========================================================
# 示例：读入两个 reconstruction（百分比 -> 比例），N=1000
# =========================================================
if __name__ == "__main__":
    # 1) 按你的文件名读取，并把百分比转为 0~1 比例
    #    你可按需替换文件名；如果只有两个，就只放两个。
    r1 = np.loadtxt('Panel_1_linking_frequency.txt')
    r2 = np.loadtxt('Panel_2_linking_frequency.txt')

    # r3 = np.loadtxt('msplMDS_linking_number.txt')
    # r4 = np.loadtxt('msplNMDS_linking_number.txt')

    percent = 100.0
    r1 = r1 / percent
    r2 = r2 / percent
    # r3 = r3 / percent
    # r4 = r4 / percent

    # 2) 名称（修正了你原来多了一个引号的小错误）
    names = ["Panel 1", "Panel 2"]
    R = [r1, r2]  # 如有更多，可扩充： R = [r1, r2, r3, r4]

    # 3) N=1000（你的设定），并设置零项阈值 const=0.005
    N = 1000
    ZERO_CONST = 0.005

    # 4) 运行：建议 strict=False，eps=max(1/(2N),1e-4) 防止 0/1 爆导数
    results = run_all_variants(
        R, names, N,
        c_range=(0.01, 10.0), n_grid=1000,
        strict=False,            # 避免 0/1 导数发散
        eps=None,                # 自动取 max(1/(2N),1e-4)
        zero_const=ZERO_CONST
    )
    pretty_print(results, zero_const=ZERO_CONST)

    # ① 看原始差值的“零项”比例（未变换）
delta_raw = r1 - r2                     # 0~1 比例空间
print("原始差 |δ|<0.005 比例:", np.mean(np.abs(delta_raw) < 0.005))

# ② 看变换后差值的“零项”比例（arcsin√p 空间）
delta_tr = np.arcsin(np.sqrt(r1)) - np.arcsin(np.sqrt(r2))
print("变换差 |δ|<0.005 比例:", np.mean(np.abs(delta_tr) < 0.005))

# ③ 看 z 的尺度（决定 u=|z|/c 的大小）
N = 1000
Var_delta_un = r1*(1-r1)/N + r2*(1-r2)/N
z_un = (delta_raw) / np.sqrt(Var_delta_un)
z_tr = (delta_tr) * np.sqrt(2*N)
print("未变换 |z| 中位数/90分位:", np.median(np.abs(z_un)), np.quantile(np.abs(z_un), 0.9))
print("变换   |z| 中位数/90分位:", np.median(np.abs(z_tr)), np.quantile(np.abs(z_tr), 0.9))
