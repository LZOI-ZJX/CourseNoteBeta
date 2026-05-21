"""
从 main.exe 生成的数据 txt 读取 u[i][j] 并绘制三维曲面。
格式：首行 n m，随后 n 行各 m 列浮点数。

用法:
  python plot.py                    # 默认 u_data.txt -> u_data.png
  python plot.py a.txt b.txt        # 每个文件生成同名 .png
"""

from pathlib import Path
import sys
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def load_u_data(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not path.is_file():
        print(f"未找到数据文件: {path}", file=sys.stderr)
        sys.exit(1)

    with path.open(encoding="utf-8") as f:
        header = f.readline().split()
        if len(header) < 2:
            sys.exit(f"首行格式错误，应为: n m（文件: {path}）")
        n, m = int(header[0]), int(header[1])
        U = np.empty((n, m), dtype=np.float64)
        for i in range(n):
            parts = f.readline().split()
            if len(parts) != m:
                sys.exit(f"第 {i + 2} 行应有 {m} 个数，实际 {len(parts)}（文件: {path}）")
            U[i, :] = [float(x) for x in parts]

    # 波动 wave_*.txt：域 (0,1)×(0,1)，结点 t_i=i/(n-1), x_j=j/(m-1)
    # 对流：t_i=i·5/n, x_j=j·10/m（与 main 对流步长一致）
    stem = path.stem.lower()
    if "wave" in stem:
        t_den = max(n - 1, 1)
        x_den = max(m - 1, 1)
        t_1d = np.arange(n, dtype=np.float64) / t_den
        x_1d = np.arange(m, dtype=np.float64) / x_den
    else:
        t_1d = np.arange(n, dtype=np.float64) * (5.0 / n)
        x_1d = np.arange(m, dtype=np.float64) * (10.0 / m)
    T, X = np.meshgrid(t_1d, x_1d, indexing="ij")
    return T, X, U


def plot_one(data_path: Path, out_path: Path | None = None) -> None:
    data_path = data_path.resolve()
    T, X, U = load_u_data(data_path)
    if out_path is None:
        out_path = data_path.with_suffix(".png")

    scheme_labels = {
        "upwind": "explicit upwind",
        "lw": "Lax–Wendroff",
        "be": "backward Euler",
    }
    stem = data_path.stem
    subtitle = ""
    for key, lab in scheme_labels.items():
        if key in stem.lower():
            subtitle = f" ({lab})"
            break

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(T, X, U, cmap="viridis", edgecolor="none", alpha=0.95)

    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$x$")
    ax.set_zlabel(r"$u$")
    ax.set_title(f"{data_path.name}{subtitle}\n{U.shape[0]}×{U.shape[1]}")

    fig.colorbar(surf, ax=ax, shrink=0.6, aspect=20, label=r"$u$")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path.resolve()}")
    plt.close(fig)


def default_data_path(script_dir: Path) -> Path:
    p = script_dir / "u_data.txt"
    if p.is_file():
        return p
    p2 = Path.cwd() / "u_data.txt"
    if p2.is_file():
        return p2
    print("未找到 u_data.txt，请传入数据文件路径。", file=sys.stderr)
    print("用法: python plot.py [file1.txt file2.txt ...]", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    paths: List[Path]
    if len(sys.argv) > 1:
        paths = [Path(a) for a in sys.argv[1:]]
    else:
        paths = [default_data_path(script_dir)]

    for p in paths:
        plot_one(p)
    if plt.get_backend().lower() != "agg" and len(paths) == 1:
        plt.show()


if __name__ == "__main__":
    main()
