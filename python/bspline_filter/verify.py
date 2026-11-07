"""One-command verification of the headline results.

    python -m bspline_filter.verify

1. Filter spectrum on a bounded grid: B-spline filter eigenvalues are real and
   in [0, 1] for every tested strength, resolution and vector strength; the
   high-order compact-filter boundary closures have amplifying, complex
   eigenvalues (paper Table 4, dissertation Table 12.4).
2. Warm-started cubic inverse for the vector filter strength: iterations per
   time step to reach the round-off floor.
3. Newton-optimized GCV: iterations to select the filter strength.
4. Filter consistency: the filtered field converges to the unfiltered one at
   high order under grid refinement (formal order is preserved).
5. Under-resolved 2-D Burgers (paper Table 9, dissertation Table 14.1):
   L-infinity error without and with the B-spline filter.
"""
from __future__ import annotations

import time

import numpy as np

from . import burgers2d
from .basis import Basis1D
from .filter import Filter1D, filter_matrix
from .gcv import ScalarGCV, SpectralDirection
from .inverse import cubic_inverse
from .stability import compact_filter_operator, penalty_pencil_eigenvalues, spectrum_summary


def check_spectrum() -> float:
    """B-spline filter spectrum over strengths and resolutions; compact closures for contrast."""
    worst_rel, worst_imag, worst_max = 0.0, 0.0, 0.0
    rng = np.random.default_rng(1)
    for n in (32, 64, 128):
        basis = Basis1D(0.0, 2.0 * np.pi, n, 7)
        B, P = basis.matrix(0), basis.scaled_penalty(6)
        strengths = [np.sqrt(x) for x in (1e-8, 1e-6, 1e-4, 1e-2, 1e-1)]
        strengths.append(0.3 * rng.random(n))
        for lam in strengths:
            nu = penalty_pencil_eigenvalues(B, P, lam)
            worst_rel = min(worst_rel, nu.min() / nu.max())
            ev = np.linalg.eigvals(filter_matrix(B, P, lam))
            worst_imag = max(worst_imag, np.abs(ev.imag).max())
            worst_max = max(worst_max, np.abs(ev).max())
    print("[1] Filter spectrum on a bounded grid (degree 7, 6th-derivative penalty)")
    print("    B-spline filter, N = 32/64/128, 5 scalar + 1 vector strength each:")
    print(f"      eigenvalues 1/(1+nu), min(nu)/max(nu) = {worst_rel:.1e} (nu >= 0 to round-off)")
    print(f"      explicit operator: max |eig| = {worst_max:.6f}, max |imag| = {worst_imag:.1e}")
    for bnd in ("8th", "6th"):
        s = spectrum_summary(compact_filter_operator(64, bnd))
        print(f"    Compact filter, {bnd}-order boundary closure: max |eig| = {s['max_mag']:.6f}, "
              f"max |imag| = {s['max_imag']:.6f}")
    return worst_rel


def check_inverse(steps: int = 100) -> dict:
    """Iterations of the warm-started cubic inverse over a slowly varying vector strength."""
    rng = np.random.default_rng(0)
    n = 80
    basis = Basis1D(0.0, 2.0 * np.pi, n, 7)
    B, P = basis.matrix(0), basis.penalty(2) * (2.0 * np.pi / n) ** 2
    lam = 1e-2 * (1.0 + 0.5 * rng.random(n))
    Q = B.T @ B + np.outer(lam, lam) * P
    V = np.linalg.inv(Q)
    floor = np.linalg.norm(np.eye(n) - Q @ V, np.inf)
    its, final = [], []
    for _ in range(steps):
        lam = lam * (1.0 + 1e-3 * rng.standard_normal(n))
        Q = B.T @ B + np.outer(lam, lam) * P
        V, hist = cubic_inverse(Q, V)
        its.append(len(hist) - 1)
        final.append(hist[-1])
    print("[2] Warm-started cubic inverse, N = 80, lambda lambda^T ~ 1e-4, 0.1% change per step")
    print(f"    iterations per step: mean {np.mean(its):.2f}, max {max(its)}; "
          f"final residual <= {max(final):.1e} (direct inverse: {floor:.1e})")
    return {"mean": float(np.mean(its)), "max": int(max(its)), "final": float(max(final))}


def check_gcv() -> int:
    """Newton iterations for the per-direction scalar GCV on an aliased 2-D field."""
    basis = Basis1D(0.0, 2.0 * np.pi, 20, 7)
    B, P = basis.matrix(0), basis.penalty(2) * (2.0 * np.pi / 20) ** 2
    X, Y = np.meshgrid(basis.x, basis.x, indexing="ij")
    rng = np.random.default_rng(2)
    u = burgers2d.exact("u", 1.0 / 80.0, X, Y, 0.5) + 1e-3 * rng.standard_normal(X.shape)
    gcv = ScalarGCV([SpectralDirection.from_operators(B, P)] * 2)
    mu, it = gcv.optimize(u)
    print("[3] Newton-optimized GCV, 20 x 20 field, one strength per direction")
    print(f"    converged in {it} Newton iterations, lambda^2 = ({mu[0]:.2e}, {mu[1]:.2e})")
    return it


def check_order() -> float:
    """Filter consistency: ||A u - u||_inf for a smooth field under grid refinement."""
    ns = [24, 32, 48, 64, 96]
    errs = []
    for n in ns:
        basis = Basis1D(0.0, 2.0 * np.pi, n, 7)
        u = np.exp(np.sin(basis.x))
        f = Filter1D(basis.matrix(0), basis.scaled_penalty(6), np.sqrt(1e-6))
        errs.append(np.max(np.abs(f.apply(u) - u)))
    rates = -np.diff(np.log(errs)) / np.diff(np.log(ns))
    print("[4] Filter consistency on a smooth field, u = exp(sin x), lambda^2 = 1e-6")
    print("    N:     " + "  ".join(f"{n:>8d}" for n in ns))
    print("    error: " + "  ".join(f"{e:8.1e}" for e in errs))
    print("    observed order: " + ", ".join(f"{r:.1f}" for r in rates))
    return float(rates.min())


def check_burgers() -> dict:
    """Table 14.1 errors for the under-resolved 2-D Burgers problem."""
    print("[5] 2-D Burgers, degree-7 B-splines, 20 x 20 points, nu = 1/80, t = 1")
    out = {}
    for case in ("unfiltered", "optimized-scalar", "optimized-vector"):
        r = burgers2d.solve(case)
        out[case] = r["linf_u"]
        print(f"    {case:<17s} L_inf(u) = L_inf(v) = {r['linf_u']:.2e}")
    print(f"    error reduction with the vector filter: "
          f"{out['unfiltered'] / out['optimized-vector']:.1f}x")
    return out


def main() -> None:
    """Run all checks and print the headline numbers."""
    t0 = time.perf_counter()
    check_spectrum()
    check_inverse()
    check_gcv()
    check_order()
    check_burgers()
    print(f"done in {time.perf_counter() - t0:.0f} s")


if __name__ == "__main__":
    main()
