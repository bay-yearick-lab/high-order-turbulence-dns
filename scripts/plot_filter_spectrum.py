"""Filter eigenvalues on a bounded grid: B-spline filter versus compact closures.

    python scripts/plot_filter_spectrum.py      # writes scripts/out/filter_spectrum.png

Reproduces the comparison behind paper Table 4: the B-spline filter spectrum is
real and inside [0, 1]; the 8th-order compact filter with 8th- and 6th-order
one-sided boundary closures has eigenvalues outside the unit disk.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from bspline_filter.stability import bspline_filter_operator, compact_filter_operator  # noqa: E402


def main():
    """Plot the eigenvalues in the complex plane."""
    plt.style.use(str(ROOT / "tools" / "lab.mplstyle"))
    n = 64
    fig, ax = plt.subplots(figsize=(7, 5))
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="#8D99AE", lw=1, label="unit circle")
    for F, lab, c, m in [(compact_filter_operator(n, "8th"), "compact, 8th-order closure", "#F2A541", "s"),
                         (compact_filter_operator(n, "6th"), "compact, 6th-order closure", "#E07A5F", "^"),
                         (bspline_filter_operator(n, lam2=1e-6), "B-spline, $\\lambda^2 = 10^{-6}$", "#0F9D8A", "o")]:
        ev = np.linalg.eigvals(F)
        ax.scatter(ev.real, ev.imag, s=22, marker=m, color=c, label=lab, zorder=3)
    ax.set_xlim(-0.1, 1.15)
    ax.set_ylim(-0.25, 0.25)
    ax.set_xlabel("Re(eigenvalue)")
    ax.set_ylabel("Im(eigenvalue)")
    ax.legend(loc="upper right")
    ax.set_title("Filter eigenvalues, 64 points, bounded grid", loc="left")
    out = ROOT / "scripts" / "out"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "filter_spectrum.png")
    print(out / "filter_spectrum.png")


if __name__ == "__main__":
    main()
