import numpy as np
from itertools import combinations

# =========================
# 基础工具（严格：默认不裁剪到 [0,1] 边界）
# =========================
def _check_open01(x, name="p"):
    if np.any((x <= 0) | (x >= 1)):
        raise ValueError(f"{name} contains values outside the open interval (0,1). "
                         f"Paper assumes 0<p<1. Found min={x.min():.4g}, max={x.max():.4g}")

def maybe_clip_open01(x, strict=True, tol=1e-12):
    if strict:
        _check_open01(x)
        return x
    # 非严格：仅在需要时极轻微裁剪到开区间（不建议；为实务稳健性）
    x = np.asarray(x, float)
    x = np.where(x <= 0, tol, x)
    x = np.where(x >= 1, 1 - tol, x)
    return x

def arcsin_sqrt(p):   # 仅在 0<p<1 定义良好
    return np.arcsin(np.sqrt(p))

def d_arcsin_sqrt_dp(p):  # d/dp arcsin(sqrt(p)) = 1 / [2 sqrt(p(1-p))]
    return 1.0 / (2.0 * np.sqrt(p * (1.0 - p)))

# =========================
# 两个平滑函数及导数（对 u 的导数）
# =========================
def G1(u):
    return (u*u) / (1.0 + u*u)
def dG1_du(u):
    return 2.0 * u / (1.0 + u*u)**2

def G2(u):
    out = np.zeros_like(u)
    m = u > 0.5
    v = (u[m] - 0.5)
    out[m] = (v*v) / (1.0 + v*v)
    return out
def dG2_du(u):
    out = np.zeros_like(u)
    m = u > 0.5
    v = (u[m] - 0.5)
    out[m] = 2.0 * v / (1.0 + v*v)**2
    return out

G_FUNCS = {"G1": (G1, dG1_du), "G2": (G2, dG2_du)}

# =========================
# δ̂、z 及其对 p、q 的导数（两种统计变换）
# =========================
def stats_untransformed(p, q, N):
    # δ̂、Var、z 以及 ∂z/∂p, ∂z/∂q
    delta = p - q
    Var = p*(1.0 - p)/N + q*(1.0 - q)/N          # Var(δ̂)
    S = np.sqrt(Var)
    z = delta / S

    # ∂Var/∂p = (1-2p)/N, ∂Var/∂q = (1-2q)/N
    dVar_dp = (1.0 - 2.0*p) / N
    dVar_dq = (1.0 - 2.0*q) / N

    # ∂z/∂p = (1/S) - δ̂ * (1/(2 S^3)) * ∂Var/∂p
    dz_dp = (1.0 / S) - delta * (0.5 * dVar_dp) / (S**3)
    # ∂z/∂q = (-1/S) - δ̂ * (1/(2 S^3)) * ∂Var/∂q
    dz_dq = (-1.0 / S) - delta * (0.5 * dVar_dq) / (S**3)

    # ∂δ̂/∂p = 1, ∂δ̂/∂q = -1
    ddelta_dp = np.ones_like(p)
    ddelta_dq = -np.ones_like(q)

    # Var(p̂)=p(1-p)/N, Var(q̂)=q(1-q)/N （供 F̂(c) term3 使用）
    Var_p = p*(1.0 - p)/N
    Var_q = q*(1.0 - q)/N

    return delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, Var_p, Var_q

def stats_transformed(p, q, N):
    tp = arcsin_sqrt(p)
    tq = arcsin_sqrt(q)
    delta = tp - tq
    # Var(delta) = 1/(2N) ⇒ z = δ̂ / sqrt(1/(2N)) = δ̂ * sqrt(2N)
    z = delta * np.sqrt(2.0 * N)

    # ∂δ̂/∂p = d_arcsin_sqrt_dp(p), ∂δ̂/∂q = -d_arcsin_sqrt_dp(q)
    ddelta_dp = d_arcsin_sqrt_dp(p)
    ddelta_dq = -d_arcsin_sqrt_dp(q)

    # z = δ̂ * sqrt(2N) ⇒ ∂z = sqrt(2N) * ∂δ̂
    dz_dp = np.sqrt(2.0 * N) * ddelta_dp
    dz_dq = np.sqrt(2.0 * N) * ddelta_dq

    # 变换情形下，F̂(c) 的校正项不再用 Var(p̂)、Var(q̂)，而是常数 1/(4N)
    return delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq

# =========================
# δ̃(c) 及其对 p、q 的解析偏导
# δ̃ = δ̂ * G(u),  u = |z|/c
# ∂δ̃/∂p = (∂δ̂/∂p) * G(u) + δ̂ * G'(u) * (1/c) * sgn(z) * (∂z/∂p)
# q 同理
# =========================
def delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du):
    u = np.abs(z) / c
    G_u = G(u)
    dG = dG_du(u)

    sign_z = np.sign(z)
    # 当 z==0 时，sgn(z)=0 合理（论文没特别指明，此处采用连续延拓）
    d_absz_dp = sign_z * dz_dp
    d_absz_dq = sign_z * dz_dq

    delta_tilde = delta * G_u

    ddt_dp = ddelta_dp * G_u + delta * dG * (1.0/c) * d_absz_dp
    ddt_dq = ddelta_dq * G_u + delta * dG * (1.0/c) * d_absz_dq
    return delta_tilde, ddt_dp, ddt_dq

# =========================
# 论文的 F̂(c)
# =========================
def Fhat_untransformed(p, q, N, c, G, dG_du, strict=True):
    p = maybe_clip_open01(np.asarray(p, float), strict=strict)
    q = maybe_clip_open01(np.asarray(q, float), strict=strict)

    delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, Var_p, Var_q = stats_untransformed(p, q, N)
    delta_tilde, ddt_dp, ddt_dq = delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du)

    term1 = np.sum(delta_tilde**2)
    term2 = -2.0 * np.sum(delta * delta_tilde)
    term3 = 2.0 * np.sum(Var_p * ddt_dp - Var_q * ddt_dq)

    return float(term1 + term2 + term3), delta_tilde

def Fhat_transformed(p, q, N, c, G, dG_du, strict=True):
    p = maybe_clip_open01(np.asarray(p, float), strict=strict)
    q = maybe_clip_open01(np.asarray(q, float), strict=strict)

    delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq = stats_transformed(p, q, N)
    delta_tilde, ddt_dp, ddt_dq = delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du)

    term1 = np.sum(delta_tilde**2)
    term2 = -2.0 * np.sum(delta * delta_tilde)
    # 变换情形的修正项： (2/(4N)) Σ{ ∂δ̃/∂p̂ - ∂δ̃/∂q̂ } = (1/(2N)) Σ{...}
    term3 = (1.0 / (2.0 * N)) * np.sum(ddt_dp - ddt_dq)

    return float(term1 + term2 + term3), delta_tilde

# =========================
# 在网格上最小化 F̂(c) → ĉ
# =========================
def _is_un(mode):
    return mode in ("un", "untransformed")

def _is_tr(mode):
    return mode in ("tr", "transformed")

def minimize_over_c(p, q, N, mode, G_name, c_min=0.01, c_max=5.0, n_grid=1000, strict=True):
    G, dG = G_FUNCS[G_name]
    cs = np.linspace(c_min, c_max, n_grid)

    best_c, best_F, best_delta_t = None, np.inf, None

    for c in cs:
        if _is_un(mode):
            F, dt = Fhat_untransformed(p, q, N, c, G, dG, strict=strict)
        elif _is_tr(mode):
            F, dt = Fhat_transformed(p, q, N, c, G, dG, strict=strict)
        else:
            raise ValueError(f"Unknown mode: {mode}")
        if F < best_F:
            best_c, best_F, best_delta_t = c, F, dt

    # 原始 delta（按当前坐标系）
    if _is_un(mode):
        delta, *_ = stats_untransformed(p, q, N)
    else:
        delta, *_ = stats_transformed(p, q, N)

    const = 0.005
    zero_ratio_tilde = float(np.mean(np.abs(best_delta_t) < const))

    # 两种基线零比例（便于和论文对齐）
    delta_un = p - q
    delta_tr = arcsin_sqrt(p) - arcsin_sqrt(q)
    baseline_un = float(np.mean(np.abs(delta_un) < const))
    baseline_tr = float(np.mean(np.abs(delta_tr) < const))

    # 诊断：u 与 G 分布
    if _is_un(mode):
        delta_mode, z_mode, *_ = stats_untransformed(p, q, N)
    else:
        delta_mode, z_mode, *_ = stats_transformed(p, q, N)
    u = np.abs(z_mode) / best_c
    G_u = G(u)
    diag = {
        "u_median": float(np.median(u)),
        "u_q90": float(np.quantile(u, 0.9)),
        "G_median": float(np.median(G_u)),
        "G_q90": float(np.quantile(G_u, 0.9)),
        "num_new_zeros": int(np.sum((np.abs(best_delta_t) < const) & (np.abs(delta_mode) >= const))),
        "num_total": int(len(delta_mode)),
    }

    return float(best_c), float(best_F), zero_ratio_tilde, baseline_un, baseline_tr, diag




# =========================
# 主流程：四种组合（严格按论文）
# =========================
def run_all_variants(R, names, N, c_range=(0.01, 10.0), n_grid=1000, strict=True):
    results = {}
    for (a, b) in combinations(range(len(R)), 2):
        p = R[a]; q = R[b]
        key = f"{names[a]} vs {names[b]}"
        results[key] = {}
        for mode in ["un", "tr"]:
            for G_name in ["G1", "G2"]:
                c_hat, Fmin, zr_tilde, base_un, base_tr, diag = minimize_over_c(
                    p, q, N, mode, G_name,
                    c_min=c_range[0], c_max=c_range[1], n_grid=n_grid, strict=strict
                )
                results[key][f"{mode}-{G_name}"] = {
                    "c_hat": c_hat,
                    "Fhat_min": Fmin,
                    "zero_ratio_tilde": zr_tilde,
                    "baseline_un": base_un,
                    "baseline_tr": base_tr,
                    "diag": diag
                }
    return results
def pretty_print(results):
    for pair, d in results.items():
        print("==== (Untransformed) ====")
        # G1
        v = d["un-G1"]
        print(f"G1: ĉ = {v['c_hat']:.4f},  F̂min = {v['Fhat_min']:.4f},  "
              f"Zero ratio (δ̃, δ) = {v['zero_ratio_tilde']*100:.2f}% , {v['baseline_un']*100:.2f}%")
        # G2
        v = d["un-G2"]
        print(f"G2: ĉ = {v['c_hat']:.4f},  F̂min = {v['Fhat_min']:.4f},  "
              f"Zero ratio (δ̃, δ) = {v['zero_ratio_tilde']*100:.2f}% , {v['baseline_un']*100:.2f}%\n")

        print("==== (Transformed) ====")
        # G1
        v = d["tr-G1"]
        print(f"G1: ĉ = {v['c_hat']:.4f},  F̂min = {v['Fhat_min']:.4f},  "
              f"Zero ratio (δ̃) = {v['zero_ratio_tilde']*100:.2f}%")
        # G2
        v = d["tr-G2"]
        print(f"G2: ĉ = {v['c_hat']:.4f},  F̂min = {v['Fhat_min']:.4f},  "
              f"Zero ratio (δ̃) = {v['zero_ratio_tilde']*100:.2f}%")



# =========================
# 读取你的 4 个 reconstruction
# =========================
r1 = np.loadtxt('panel_1_max.txt')
r2 = np.loadtxt('panel_2_max.txt')


names = ["Panel 1", "Panel 2"]

percent = 100.0
r1 = r1/ percent
r2 = r2/ percent

R = [r1, r2]

# 注意：严格模式下，若存在 0 或 1，会报错；若你确实有 0/1，请先把 strict=False 运行或在读取前做极轻微处理
N = 1000
results = run_all_variants(R, names, N, c_range=(-10.0, 10.0), n_grid=10000, strict=True)
pretty_print(results)

import matplotlib.pyplot as plt
import os

def scan_and_save_F(p, q, N, mode, G_name, save_dir=".",
                    c_min=-10.0, c_max=10.0, n_grid=1000, strict=True):
    G, dG = G_FUNCS[G_name]
    cs = np.linspace(c_min, c_max, n_grid)
    Fs = []
    for c in cs:
        if mode in ("un", "untransformed"):
            F, _ = Fhat_untransformed(p, q, N, c, G, dG, strict=strict)
        else:
            F, _ = Fhat_transformed(p, q, N, c, G, dG, strict=strict)
        Fs.append(F)
    Fs = np.asarray(Fs)

    best_idx = np.argmin(Fs)
    best_c = float(cs[best_idx])
    best_F = float(Fs[best_idx])
    print(f"[{mode}-{G_name}] best c ≈ {best_c:.4f},  Fmin = {best_F:.6f}")

    import matplotlib.pyplot as plt
    import os, time
    plt.figure(figsize=(7,5))
    plt.plot(cs, Fs, '-', lw=1.5)
    plt.axvline(best_c, color='r', ls='--', label=f"ĉ = {best_c:.4f}")
    plt.xlabel("c (smoothing parameter)")
    plt.ylabel("F̂(c)")
    plt.title(f"F̂(c) curve: {mode}-{G_name}")
    plt.grid(True, alpha=0.3)
    plt.legend()

    os.makedirs(save_dir, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    # 在文件名里带上 best_c，避免看错旧图
    filename = f"Fhat_{mode}_{G_name}_c{best_c:.4f}_{ts}.png"
    save_path = os.path.join(save_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ 图像已保存到 {save_path}")


# 保存到当前目录 /app 下
out_dir = "/app"

scan_and_save_F(r1, r2, N, "untransformed", "G1", save_dir=out_dir)
scan_and_save_F(r1, r2, N, "untransformed", "G2", save_dir=out_dir)
scan_and_save_F(r1, r2, N, "transformed", "G1", save_dir=out_dir)
scan_and_save_F(r1, r2, N, "transformed", "G2", save_dir=out_dir)

