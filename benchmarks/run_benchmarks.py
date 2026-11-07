"""Timing and iteration benchmarks for the B-spline filter components.

    python benchmarks/run_benchmarks.py            # prints a Markdown table

1. 3-D filter application on a 60 x 80 x 60 channel-size field with Kronecker
   line operators (banded LAPACK solves), compared with the coefficient solve
   B^{-1} u that every collocation right-hand side already performs.
2. Newton-optimized GCV, one strength per direction, on the same field size.
3. Warm-started cubic inverse for a slowly varying vector strength.
"""
from __future__ import annotations

import platform
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from bspline_filter.basis import Basis1D  # noqa: E402
from bspline_filter.filter import BandedLU, Filter1D, filter_nd  # noqa: E402
from bspline_filter.gcv import ScalarGCV, SpectralDirection  # noqa: E402
from bspline_filter.inverse import cubic_inverse  # noqa: E402


def best_of(fn, repeat=5):
    """Minimum wall time of ``repeat`` calls."""
    ts = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return min(ts)


def main() -> None:
    """Run the benchmarks and print a Markdown table."""
    shape = (60, 80, 60)
    bases = [Basis1D(0.0, L, n, 7) for n, L in zip(shape, (4 * np.pi, 2.0, 4 * np.pi / 3))]
    Bs = [b.matrix(0) for b in bases]
    Ps = [b.scaled_penalty(6) for b in bases]
    rng = np.random.default_rng(0)
    u = rng.standard_normal(shape)
    rows = []

    lus = [BandedLU(B) for B in Bs]

    def coeff_solve():
        v = u
        for ax, lu in enumerate(lus):
            w = np.moveaxis(v, ax, 0)
            s = w.shape
            v = np.moveaxis(lu.solve(np.asfortranarray(w.reshape(s[0], -1))).reshape(s), 0, ax)
        return v

    filters = [Filter1D(B, P, 1e-3) for B, P in zip(Bs, Ps)]
    t_coef = best_of(coeff_solve)
    t_filt = best_of(lambda: filter_nd(u, filters))
    rows.append(("Coefficient solve B^-1 u, 60x80x60 (3 banded line sweeps)", f"{t_coef * 1e3:.1f} ms", ""))
    rows.append(("B-spline filter A u, 60x80x60 (3 banded line sweeps)", f"{t_filt * 1e3:.1f} ms",
                 f"{t_filt / t_coef:.2f}x the coefficient solve"))

    dirs = [SpectralDirection.from_operators(B, P) for B, P in zip(Bs, Ps)]
    gcv = ScalarGCV(dirs)
    noisy = np.sin(bases[0].x)[:, None, None] * np.cos(np.pi * bases[1].x)[None, :, None] \
        * np.ones(shape) + 1e-3 * rng.standard_normal(shape)
    mu, it = gcv.optimize(noisy)
    rows.append(("Newton-optimized GCV, 3 strengths, 60x80x60", f"{it} Newton iterations",
                 f"lambda^2 = {mu[0]:.1e}, {mu[1]:.1e}, {mu[2]:.1e}"))

    n = 80
    b = Basis1D(0.0, 2.0 * np.pi, n, 7)
    B, P = b.matrix(0), b.penalty(2) * (2.0 * np.pi / n) ** 2
    lam = 1e-2 * (1.0 + 0.5 * rng.random(n))
    V = np.linalg.inv(B.T @ B + np.outer(lam, lam) * P)
    its, t_it = [], 0.0
    for _ in range(200):
        lam = lam * (1.0 + 1e-3 * rng.standard_normal(n))
        Q = B.T @ B + np.outer(lam, lam) * P
        t0 = time.perf_counter()
        V, h = cubic_inverse(Q, V)
        t_it += time.perf_counter() - t0
        its.append(len(h) - 1)
    rows.append(("Warm-started cubic inverse, N = 80, 200 steps",
                 f"{t_it / 200 * 1e3:.2f} ms/step", f"mean {np.mean(its):.2f}, max {max(its)} iterations"))

    print(f"Hardware: {platform.processor() or platform.machine()}, Python {platform.python_version()}, "
          f"NumPy {np.__version__}\n")
    print("| Benchmark | Time | Notes |")
    print("| --- | --- | --- |")
    for r in rows:
        print(f"| {r[0]} | {r[1]} | {r[2]} |")


if __name__ == "__main__":
    main()
