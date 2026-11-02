import numpy as np
from itertools import combinations
from pathlib import Path

# =========================================================
# Utilities: validate boundaries and optionally clip values (the paper requires 0 < p < 1, or derivatives diverge).
# =========================================================
def _check_open01(x, name="p"):
    if np.any((x <= 0) | (x >= 1)):
        raise ValueError(
            f"{name} 含有不在开区间(0,1)内的值。论文假设 0<p<1。"
            f" 发现 min={x.min():.4g}, max={x.max():.4g}"
        )

def maybe_clip_open01(x, strict=True, eps=1e-4):
    """
    strict=True: Require 0 < p < 1 and raise an error at the boundaries (the paper's original assumption).
    strict=False: For numerical robustness, minimally clip values to (eps, 1-eps).
    Recommended: strict=False with eps at least 1/(2N).
    """
    x = np.asarray(x, float)
    if strict:
        _check_open01(x)
        return x
    # Non-strict mode: apply minimal clipping.
    low = eps
    high = 1.0 - eps
    return np.clip(x, low, high)

# =========================================================
# Transformation and derivative (arcsin square-root transform), used in the paper for variance stabilization.
# =========================================================
def arcsin_sqrt(p):
    # Well-defined only for 0 < p < 1.
    return np.arcsin(np.sqrt(p))

def d_arcsin_sqrt_dp(p):
    # d/dp arcsin(sqrt(p)) = 1 / [2 sqrt(p(1-p))]
    return 1.0 / (2.0 * np.sqrt(p * (1.0 - p)))

# =========================================================
# Two G functions for semi-soft thresholding and their derivatives with respect to u.
# u = |z|/c
# =========================================================
def G1(u):
    # Smooth shrinkage: strong for small u and approaching 1 for large u.
    return (u*u) / (1.0 + u*u)

def dG1_du(u):
    return 2.0 * u / (1.0 + u*u)**2

def G2(u):
    # Includes a hard threshold: 0 for u <= 0.5, then a smooth function for u > 0.5.
    out = np.zeros_like(u)
    m = u > 0.5
    v = (u[m] - 0.5)
    out[m] = (v*v) / (1.0 + v*v)
    return out

def dG2_du(u):
    # The derivative is 0 for u <= 0.5 (zero is a valid choice of subgradient).
    out = np.zeros_like(u)
    m = u > 0.5
    v = (u[m] - 0.5)
    out[m] = 2.0 * v / (1.0 + v*v)**2
    return out

G_FUNCS = {"G1": (G1, dG1_du), "G2": (G2, dG2_du)}

# =========================================================
# Statistics for the untransformed and transformed approaches.
# Returns:
#   delta: difference (p-q when untransformed; arcsin(sqrt(p))-arcsin(sqrt(q)) when transformed)
#   z: standardized difference based on the corresponding variance approximation
#   ddelta_dp, ddelta_dq: derivatives of delta with respect to p and q
#   dz_dp, dz_dq: derivatives of z with respect to p and q (used for the Stein correction)
#   Var_p, Var_q when using the untransformed approach
# Assumption: as requested, Np=Nq=N=1000.
# =========================================================
def stats_untransformed(p, q, N):
    """
    Untransformed approach:
        delta = p - q
        Var(delta) = p(1-p)/N + q(1-q)/N
        z = delta / sqrt(Var(delta))
    """
    delta = p - q
    Var = p*(1.0 - p)/N + q*(1.0 - q)/N
    S = np.sqrt(Var)
    z = delta / S

    # Derivatives of Var with respect to p and q.
    dVar_dp = (1.0 - 2.0*p) / N
    dVar_dq = (1.0 - 2.0*q) / N

    # Derivatives of z with respect to p and q, where z = delta / S.
    dz_dp = (1.0 / S) - delta * (0.5 * dVar_dp) / (S**3)
    dz_dq = (-1.0 / S) - delta * (0.5 * dVar_dq) / (S**3)

    # Derivatives of delta with respect to p and q.
    ddelta_dp = np.ones_like(p)
    ddelta_dq = -np.ones_like(q)

    # Values used by the Stein correction.
    Var_p = p*(1.0 - p)/N
    Var_q = q*(1.0 - q)/N

    return delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, Var_p, Var_q

def stats_transformed(p, q, N):
    """
    Transformed (arcsin(sqrt(p))) approach:
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
# Shrunk difference delta_tilde and its derivatives.
#   delta_tilde = delta * G(u), where u = |z|/c
#   d(delta_tilde)/dp = (d(delta)/dp)G(u) + delta * G'(u) * (1/c) * sign(z) * dz/dp
#   The derivative with respect to q is analogous.
# =========================================================
def delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du):
    u = np.abs(z) / c
    G_u = G(u)
    dG = dG_du(u)

    sign_z = np.sign(z)  # Use 0 at z=0, consistent with the continuous extension.
    d_absz_dp = sign_z * dz_dp
    d_absz_dq = sign_z * dz_dq

    delta_tilde = delta * G_u

    ddt_dp = ddelta_dp * G_u + delta * dG * (1.0/c) * d_absz_dp
    ddt_dq = ddelta_dq * G_u + delta * dG * (1.0/c) * d_absz_dq
    return delta_tilde, ddt_dp, ddt_dq

# =========================================================
# The paper's F_hat(c) criterion for both approaches.
# =========================================================
def Fhat_untransformed(p, q, N, c, G, dG_du, strict=True, eps=1e-4):
    # Convert inputs to proportions and handle 0/1 according to strict and eps.
    p = maybe_clip_open01(np.asarray(p, float), strict=strict, eps=eps)
    q = maybe_clip_open01(np.asarray(q, float), strict=strict, eps=eps)

    delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, Var_p, Var_q = stats_untransformed(p, q, N)
    delta_tilde, ddt_dp, ddt_dq = delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du)

    term1 = np.sum(delta_tilde**2)
    term2 = -2.0 * np.sum(delta * delta_tilde)
    # Stein correction for the untransformed approach.
    term3 = 2.0 * np.sum(Var_p * ddt_dp - Var_q * ddt_dq)

    return float(term1 + term2 + term3), delta_tilde

def Fhat_transformed(p, q, N, c, G, dG_du, strict=True, eps=1e-4):
    p = maybe_clip_open01(np.asarray(p, float), strict=strict, eps=eps)
    q = maybe_clip_open01(np.asarray(q, float), strict=strict, eps=eps)

    delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq = stats_transformed(p, q, N)
    delta_tilde, ddt_dp, ddt_dq = delta_tilde_and_partials(delta, z, ddelta_dp, ddelta_dq, dz_dp, dz_dq, c, G, dG_du)

    term1 = np.sum(delta_tilde**2)
    term2 = -2.0 * np.sum(delta * delta_tilde)
    # Transformed-case Stein correction: (2/(4N)) * sum(d(delta_tilde)/dp_hat - d(delta_tilde)/dq_hat) = (1/(2N)) * sum(...).
    term3 = (1.0 / (2.0 * N)) * np.sum(ddt_dp - ddt_dq)

    return float(term1 + term2 + term3), delta_tilde

# =========================================================
# Minimize F_hat(c) over a grid to obtain c_hat; calculate zero_ratio using the paper's 0.005 threshold.
# =========================================================
def minimize_over_c(p, q, N, mode, G_name,
                    c_min=0.001, c_max=60.0, n_grid=1000,
                    strict=False, eps=1e-4, zero_const=0.005):
    """
    Parameters:
      mode: "untransformed" or "transformed"
      G_name: "G1" or "G2"
      zero_const: threshold recommended by the paper (default 0.005)
      strict=False: minimally clips values by default to prevent divergent derivatives at 0/1
      eps: clipping amount in non-strict mode (recommended to be >= 1/(2N))
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

    # The paper's zero-term criterion: |delta_tilde| < const (const=0.005).
    zero_ratio = float(np.mean(np.abs(best_delta_t) < zero_const))
    return float(best_c), float(best_F), zero_ratio

# =========================================================
# Main workflow: compare each pair of reconstructions and report all four combinations.
# =========================================================
def run_all_variants(R, names, N,
                     c_range=(0.01, 10.0), n_grid=1000,
                     strict=False, eps=None, zero_const=0.005):
    """
    R: list of proportion vectors (0-1), one for each reconstruction
    names: list of labels, one for each vector
    N: sample size (assumed constant at N_i = N = 1000)
    eps: clipping amount in non-strict mode; defaults to max(1/(2N), 1e-4)
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
# Example: load two reconstructions (percentages to proportions), with N=1000.
# =========================================================
if __name__ == "__main__":
    # 1) Load the named files and convert percentages to proportions from 0 to 1.
    #    Replace the filenames as needed; include only the available inputs.
    script_dir = Path(__file__).resolve().parent
    r1 = np.loadtxt(script_dir / 'panel_1_max.txt')
    r2 = np.loadtxt(script_dir / 'panel_2_max.txt')

    # r3 = np.loadtxt('msplMDS_linking_number.txt')
    # r4 = np.loadtxt('msplNMDS_linking_number.txt')

    percent = 100.0
    r1 = r1 / percent
    r2 = r2 / percent
    # r3 = r3 / percent
    # r4 = r4 / percent

    # 2) Labels.
    names = ["Panel 1", "Panel 2"]
    R = [r1, r2]  # Add more reconstructions here if needed: R = [r1, r2, r3, r4].

    # 3) Set N=1000 and the zero-term threshold to 0.005.
    N = 1000
    ZERO_CONST = 0.005

    # 4) Run the analysis; strict=False and eps=max(1/(2N), 1e-4) avoid divergent derivatives at 0/1.
    results = run_all_variants(
        R, names, N,
        c_range=(0.01, 10.0), n_grid=1000,
        strict=False,            # Avoid divergent derivatives at 0/1.
        eps=None,                # Automatically use max(1/(2N), 1e-4).
        zero_const=ZERO_CONST
    )
    pretty_print(results, zero_const=ZERO_CONST)

    # 1) Measure the proportion of near-zero raw differences (untransformed).
delta_raw = r1 - r2                     # Proportion scale from 0 to 1.
print("原始差 |δ|<0.005 比例:", np.mean(np.abs(delta_raw) < 0.005))

# 2) Measure the proportion of near-zero differences after the arcsin(sqrt(p)) transform.
delta_tr = np.arcsin(np.sqrt(r1)) - np.arcsin(np.sqrt(r2))
print("变换差 |δ|<0.005 比例:", np.mean(np.abs(delta_tr) < 0.005))

# 3) Inspect the scale of z, which determines u=|z|/c.
N = 1000
Var_delta_un = r1*(1-r1)/N + r2*(1-r2)/N
z_un = (delta_raw) / np.sqrt(Var_delta_un)
z_tr = (delta_tr) * np.sqrt(2*N)
print("未变换 |z| 中位数/90分位:", np.median(np.abs(z_un)), np.quantile(np.abs(z_un), 0.9))
print("变换   |z| 中位数/90分位:", np.median(np.abs(z_tr)), np.quantile(np.abs(z_tr), 0.9))
