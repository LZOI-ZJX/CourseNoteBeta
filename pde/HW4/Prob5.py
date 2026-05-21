

from __future__ import annotations

import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (OSError, ValueError):
        pass

import matplotlib
import matplotlib.pyplot as plt


def initial_profile(x: np.ndarray) -> np.ndarray:
    out = np.empty_like(x, dtype=np.float64)
    n = x.size
    for i in range(n):
        xi = float(x[i])
        out[i] = 1.0 if 0.4 <= xi <= 0.6 else 0.0
    return out


def exact_solution(x: np.ndarray, t: float) -> np.ndarray:
    """u(x,t) = u_0(x - t)，a = 1。"""
    out = np.empty_like(x, dtype=np.float64)
    n = x.size
    for i in range(n):
        xi = float(x[i]) - t
        out[i] = 1.0 if 0.4 <= xi <= 0.6 else 0.0
    return out


def step_upwind(u: np.ndarray, lam: float) -> np.ndarray:
    """迎风格式（a=1）：后向差分，周期边界。"""
    nx = int(u.size)
    u_new = np.empty_like(u, dtype=np.float64)
    for j in range(nx):
        jm = j - 1 if j > 0 else nx - 1
        u_new[j] = u[j] - lam * (u[j] - u[jm])
    return u_new


def step_lax_friedrichs(u: np.ndarray, lam: float) -> np.ndarray:
    nx = int(u.size)
    u_new = np.empty_like(u, dtype=np.float64)
    for j in range(nx):
        jm = j - 1 if j > 0 else nx - 1
        jp = j + 1 if j < nx - 1 else 0
        ul, uj, ur = float(u[jm]), float(u[j]), float(u[jp])
        u_new[j] = 0.5 * (ur + ul) - 0.5 * lam * (ur - ul)
    return u_new


def step_lax_wendroff(u: np.ndarray, lam: float) -> np.ndarray:
    nx = int(u.size)
    u_new = np.empty_like(u, dtype=np.float64)
    for j in range(nx):
        jm = j - 1 if j > 0 else nx - 1
        jp = j + 1 if j < nx - 1 else 0
        ul, uj, ur = float(u[jm]), float(u[j]), float(u[jp])
        u_new[j] = (
            uj
            - 0.5 * lam * (ur - ul)
            + 0.5 * lam**2 * (ur - 2.0 * uj + ul)
        )
    return u_new


def integrate_to_time(
    u0: np.ndarray,
    lam: float,
    n_steps: int,
    step_fn,
) -> np.ndarray:
    u = u0.copy()
    for _ in range(n_steps):
        u = step_fn(u, lam)
    return u


def print_cfl_analysis(a: float, tau: float, h: float) -> None:
    cfl = abs(a) * tau / h
    print("=== CFL 条件分析 ===")
    print(f"方程：u_t + a u_x = 0，波速 a = {a:g}")
    print(f"空间步长 h = {h:g}，时间步长 τ = {tau:g}，网格比 λ = τ/h = {cfl:g}")
    print(f"CFL 数 ν = |a| τ / h = {cfl:g}")
    print(
        "显式迎风格式、Lax–Friedrichs、Lax–Wendroff 在线性常系数情形下"
        "通常要求 ν ≤ 1（Lax–Friedrichs 在 ν=1 时退化为特殊情况）。"
    )
    if cfl <= 1.0 + 1e-14:
        print(f"当前 ν = {cfl:g} ≤ 1：满足常见线性稳定性 CFL 上界。")
    else:
        print(f"当前 ν = {cfl:g} > 1：不满足 ν ≤ 1，显式格式可能不稳定。")


def main() -> None:
    a = 1.0
    h = 0.1
    lam = 0.5
    tau = lam * h
    x_min, x_max = 0.0, 2.0
    nx = int(round((x_max - x_min) / h)) + 1
    x = np.empty(nx, dtype=np.float64)
    for i in range(nx):
        x[i] = x_min + i * h
    h_actual = float(x[1] - x[0])
    if abs(h_actual - h) > 1e-12:
        raise RuntimeError(f"网格间距与预期不符：h_actual={h_actual}, h={h}")

    print_cfl_analysis(a, tau, h)

    u0 = initial_profile(x)
    t_targets = (0.5, 1.0)
    schemes = (
        ("Upwind", step_upwind),
        ("Lax–Friedrichs", step_lax_friedrichs),
        ("Lax–Wendroff", step_lax_wendroff),
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    colors = {"Upwind": "C0", "Lax–Friedrichs": "C1", "Lax–Wendroff": "C2"}

    for ax, t_end in zip(axes, t_targets):
        n_steps = int(round(t_end / tau))
        if abs(n_steps * tau - t_end) > 1e-12:
            raise RuntimeError(f"t_end={t_end} 不能由 τ 的整数步到达")
        ax.plot(x, exact_solution(x, t_end), "k-", lw=2.2, label="Exact")
        for name, step_fn in schemes:
            u_num = integrate_to_time(u0, lam, n_steps, step_fn)
            ax.plot(
                x,
                u_num,
                "--",
                color=colors[name],
                lw=1.6,
                label=name,
            )
        ax.set_xlabel("$x$")
        ax.set_ylabel("$u$")
        ax.set_title(f"$t = {t_end:g}$")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        r"$u_t + u_x = 0$: Upwind vs Lax–Friedrichs vs Lax–Wendroff "
        rf"($h={h}$, $\lambda=\tau/h={lam}$)"
    )
    fig.tight_layout()
    out_path = "advection_comparison.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\n图像已保存：{out_path}")
    if matplotlib.get_backend().lower() == "agg":
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
