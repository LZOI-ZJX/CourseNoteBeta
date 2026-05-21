"""
Numerical analysis for main.exe output (same grid as plot.py).

- Advection: exact u=t^2*sin(x), Linf/RMS error, CFL=m/(2n).
- Wave (filename contains 'wave'): BC/IC checks, CFL=dt/dx, interior norms,
  discrete energy trend, optional cross-file ||U_a-U_b|| if same grid.
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any, Dict, List, Tuple

import numpy as np


def load_grid(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    with path.open(encoding="utf-8") as f:
        header = f.readline().split()
        if len(header) < 2:
            raise ValueError(f"bad header n m: {path}")
        n, m = int(header[0]), int(header[1])
        U = np.empty((n, m), dtype=np.float64)
        for i in range(n):
            parts = f.readline().split()
            if len(parts) != m:
                raise ValueError(f"row {i + 2} expected {m} values: {path}")
            U[i, :] = [float(x) for x in parts]

    is_wave = "wave" in path.stem.lower()
    if is_wave:
        td, xd = max(n - 1, 1), max(m - 1, 1)
        t_1d = np.arange(n, dtype=np.float64) / td
        x_1d = np.arange(m, dtype=np.float64) / xd
    else:
        t_1d = np.arange(n, dtype=np.float64) * (5.0 / n)
        x_1d = np.arange(m, dtype=np.float64) * (10.0 / m)

    return U, t_1d, x_1d, is_wave


def exact_advection(t: np.ndarray, x: np.ndarray) -> np.ndarray:
    T, X = np.meshgrid(t, x, indexing="ij")
    return (T**2) * np.sin(X)


def discrete_energy_series(U: np.ndarray, dt: float, dx: float) -> np.ndarray:
    """Rough semi-discrete energy E^n ~ sum 0.5*v^2 + 0.5*u_x^2 on interior, n=1..nt-2."""
    nt, nx = U.shape
    if nt < 3 or nx < 3:
        return np.array([])
    inv_2dt = 1.0 / (2.0 * dt)
    inv_2dx = 1.0 / (2.0 * dx)
    energies = []
    for n in range(1, nt - 1):
        v = (U[n + 1, 1:-1] - U[n - 1, 1:-1]) * inv_2dt
        ux = (U[n, 2:] - U[n, 0:-2]) * inv_2dx
        e = 0.5 * np.mean(v**2) + 0.5 * np.mean(ux**2)
        energies.append(e)
    return np.asarray(energies)


def analyze_one(path: Path) -> Dict[str, Any]:
    path = path.resolve()
    U, t_1d, x_1d, is_wave = load_grid(path)
    n, m = U.shape
    out: Dict[str, Any] = {"file": str(path), "n": n, "m": m, "is_wave": is_wave}

    finite = np.isfinite(U)
    with np.errstate(over="ignore", invalid="ignore"):
        out["u_min"] = float(np.min(U[finite])) if np.any(finite) else float("nan")
        out["u_max"] = float(np.max(U[finite])) if np.any(finite) else float("nan")
        out["u_mean_abs"] = (
            float(np.mean(np.abs(U[finite]))) if np.any(finite) else float("nan")
        )
        mean_sq = np.mean(U[finite] ** 2) if np.any(finite) else float("nan")
    out["u_rms"] = float(np.sqrt(mean_sq)) if np.isfinite(mean_sq) else float("inf")
    out["has_nan"] = bool(np.any(np.isnan(U)))
    out["has_inf"] = bool(np.any(np.isinf(U)))

    if is_wave:
        out["problem"] = "wave u_tt-u_xx=f on (0,1)^2"
        nt, nx = n, m
        dt = 1.0 / max(nt - 1, 1)
        dx = 1.0 / max(nx - 1, 1)
        out["dt"] = dt
        out["dx"] = dx
        out["cfl_wave"] = float(dt / dx)
        out["bc_left_max_abs"] = float(np.max(np.abs(U[:, 0])))
        out["bc_right_max_err"] = float(
            np.max(np.abs(U[:, -1] - np.sin(t_1d)))
        )
        out["ic_u_max_abs"] = float(np.max(np.abs(U[0, :])))
        if nt >= 2 and nx > 2:
            v_num = (U[1, 1:-1] - U[0, 1:-1]) / dt
            xmid = x_1d[1:-1]
            out["ic_ut_max_err"] = float(np.max(np.abs(v_num - xmid)))
        else:
            out["ic_ut_max_err"] = float("nan")
        interior = U[:, 1:-1] if nx > 2 else U
        rms_rows = np.sqrt(np.mean(interior**2, axis=1))
        out["rms_interior_mean_over_t"] = float(np.mean(rms_rows))
        out["rms_interior_t_final"] = float(rms_rows[-1])
        es = discrete_energy_series(U, dt, dx)
        if es.size > 0:
            out["energy_mean"] = float(np.mean(es))
            out["energy_min"] = float(np.min(es))
            out["energy_max"] = float(np.max(es))
            out["energy_rel_range"] = float((np.max(es) - np.min(es)) / (np.mean(es) + 1e-30))
        else:
            out["energy_mean"] = float("nan")
            out["energy_min"] = float("nan")
            out["energy_max"] = float("nan")
            out["energy_rel_range"] = float("nan")
        return out

    out["problem"] = "advection, exact u=t^2*sin(x)"
    dt = 5.0 / n
    dx = 10.0 / m
    out["dt"] = dt
    out["dx"] = dx
    out["cfl_a1"] = float(dt / dx)

    Ue = exact_advection(t_1d, x_1d)
    err = U - Ue
    ef = np.isfinite(err)
    with np.errstate(over="ignore", invalid="ignore"):
        if np.any(ef):
            out["err_Linf"] = float(np.max(np.abs(err[ef])))
            ms_e = np.mean(err[ef] ** 2)
            out["err_rms"] = float(np.sqrt(ms_e)) if np.isfinite(ms_e) else float("inf")
        else:
            out["err_Linf"] = float("inf")
            out["err_rms"] = float("inf")

    denom = max(float(np.max(np.abs(Ue))), 1e-30)
    rms_ue = float(np.sqrt(np.mean(Ue**2)))
    out["err_rel_Linf"] = out["err_Linf"] / denom if np.isfinite(out["err_Linf"]) else float("inf")
    out["err_rel_rms"] = (
        out["err_rms"] / max(rms_ue, 1e-30) if np.isfinite(out["err_rms"]) else float("inf")
    )

    if out["cfl_a1"] > 1.0 + 1e-9:
        out["stability_note"] = (
            f"CFL={out['cfl_a1']:.4g}>1: explicit upwind/LW violate stability; "
            "backward Euler can still be accurate (see *_be.txt)."
        )
    else:
        out["stability_note"] = f"CFL={out['cfl_a1']:.4g}<=1 (explicit schemes stable)."

    return out


def print_report(d: Dict[str, Any]) -> None:
    print(f"\n=== {d['file']} ===")
    print(f"  grid: n={d['n']} (t), m={d['m']} (x)  |  {d['problem']}")
    print(
        f"  u: min={d['u_min']:.6e}, max={d['u_max']:.6e}, RMS(u)={d['u_rms']:.6e}"
    )
    if d["has_nan"] or d["has_inf"]:
        print(f"  warning: NaN={d['has_nan']}, Inf={d['has_inf']}")

    if d.get("cfl_wave") is not None:
        print(
            f"  dt={d['dt']:.6e}, dx={d['dx']:.6e}, CFL_wave=dt/dx={d['cfl_wave']:.6f}"
        )
        if d["cfl_wave"] > 1.0 + 1e-9:
            print("  note: CFL>1 may destabilize explicit central / compact leapfrog.")
        print("  BC/IC residuals (should be ~0 or discretization size):")
        print(f"    max|u(0,t)| on x=0     = {d['bc_left_max_abs']:.6e}")
        print(f"    max|u(1,t)-sin(t)|     = {d['bc_right_max_err']:.6e}")
        print(f"    max|u(x,0)| on t=0     = {d['ic_u_max_abs']:.6e}")
        print(f"    max|D_t u(x,0)-x| (1st step vs u_t(x,0)=x) = {d['ic_ut_max_err']:.6e}")
        print("  interior RMS in x, mean over t / at t=1:")
        print(f"    mean_t RMS_x(u) = {d['rms_interior_mean_over_t']:.6e}")
        print(f"    RMS_x(u) at final t = {d['rms_interior_t_final']:.6e}")
        if np.isfinite(d["energy_mean"]):
            print(
                f"  discrete energy (interior): mean={d['energy_mean']:.6e}, "
                f"min={d['energy_min']:.6e}, max={d['energy_max']:.6e}, "
                f"rel range={d['energy_rel_range']:.6e}"
            )
        return

    if d.get("dt") is not None:
        print(f"  dt=5/n={d['dt']:.6e}, dx=10/m={d['dx']:.6e}, CFL(a=1)={d['cfl_a1']:.6f}")
        print(f"  {d.get('stability_note', '')}")
        print("  error vs exact u=t^2*sin(x):")
        print(f"    ||e||_inf = {d['err_Linf']:.6e}")
        print(f"    ||e||_RMS = {d['err_rms']:.6e}")
        print(f"    rel ||e||_inf / ||u_ex||_inf ~ {d['err_rel_Linf']:.6e}")
        print(f"    rel RMS(e)/RMS(u_ex) ~ {d['err_rel_rms']:.6e}")


def pairwise_wave_diff(paths: List[Path], results: List[Dict[str, Any]]) -> None:
    wave_idx = [i for i, r in enumerate(results) if r.get("is_wave")]
    if len(wave_idx) < 2:
        return
    loaded: List[Tuple[str, np.ndarray]] = []
    for i in wave_idx:
        U, _, _, _ = load_grid(Path(results[i]["file"]))
        loaded.append((results[i]["file"], U))
    shape0 = loaded[0][1].shape
    if not all(U.shape == shape0 for _, U in loaded):
        print("\n=== wave cross-compare: skip (grid shapes differ) ===")
        return
    print("\n=== wave cross-compare (same nt x nx) ===")
    print(f"  shape {shape0[0]} x {shape0[1]}")
    for i in range(len(loaded)):
        for j in range(i + 1, len(loaded)):
            d = loaded[i][1] - loaded[j][1]
            linf = float(np.max(np.abs(d)))
            rms = float(np.sqrt(np.mean(d**2)))
            ni = Path(loaded[i][0]).name
            nj = Path(loaded[j][0]).name
            print(f"  ||{ni} - {nj}||_inf = {linf:.6e},  RMS = {rms:.6e}")


def main() -> None:
    paths: List[Path]
    if len(sys.argv) > 1:
        paths = [Path(a) for a in sys.argv[1:]]
    else:
        here = Path(__file__).resolve().parent
        cand = [
            here / "wave_m1.txt",
            here / "wave_m2.txt",
            here / "wave_m3.txt",
            here / "wave_u.txt",
            here / "T1Result" / "u_250x1000_upwind.txt",
            here / "T1Result" / "u_250x1000_lw.txt",
            here / "T1Result" / "u_250x1000_be.txt",
            here / "u_data.txt",
        ]
        paths = [p for p in cand if p.is_file()]
        if not paths:
            print(
                "usage: python numerical_analysis.py <data1.txt> [data2.txt ...]",
                file=sys.stderr,
            )
            sys.exit(1)

    results: List[Dict[str, Any]] = []
    for p in paths:
        try:
            r = analyze_one(p)
            results.append(r)
            print_report(r)
        except (OSError, ValueError) as e:
            print(f"skip {p}: {e}", file=sys.stderr)

    pairwise_wave_diff(paths, results)


if __name__ == "__main__":
    main()
