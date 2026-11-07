"""Eigenvalue (linear stability) analysis of filters on a bounded grid.

A filter F applied every step is linearly stable only if no eigenvalue has
modulus above 1, and it is free of dispersion error only if every eigenvalue is
real. The B-spline filter A = B Q^{-1} B^T is symmetric positive semi-definite
with spectrum in [0, 1] by construction. For comparison, this module builds the
8th-order tridiagonal compact (Pade-type) filter with its one-sided 8th-, 6th-
and 4th-order boundary closures and the 2nd-order centered near-wall closure
(coefficients as tabulated in the paper, Tables 1 and 3).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from .basis import Basis1D
from .filter import filter_matrix


def penalty_pencil_eigenvalues(B: np.ndarray, P: np.ndarray, lam) -> np.ndarray:
    """Eigenvalues nu of the symmetric-definite pencil Omega v = nu B^T B v.

    Omega = lambda^2 P (scalar) or diag(lambda) P diag(lambda) (vector) is
    positive semi-definite, so nu >= 0.
    """
    lam = np.asarray(lam, dtype=float)
    omega = lam**2 * P if lam.ndim == 0 else np.outer(lam, lam) * P
    return eigh(omega, B.T @ B, eigvals_only=True)


def bspline_filter_eigenvalues(B: np.ndarray, P: np.ndarray, lam) -> np.ndarray:
    """Eigenvalues of A = B Q^{-1} B^T.

    For invertible B, A = (I + B^{-T} Omega B^{-1})^{-1}, so its eigenvalues are
    1 / (1 + nu) with nu from ``penalty_pencil_eigenvalues``; all lie in (0, 1].
    """
    return 1.0 / (1.0 + penalty_pencil_eigenvalues(B, P, lam))


def spectrum_summary(F: np.ndarray) -> dict:
    """Min/max real part, imaginary part and modulus of the eigenvalues of F."""
    ev = np.linalg.eigvals(F)
    return {"min_real": ev.real.min(), "max_real": ev.real.max(),
            "min_imag": ev.imag.min(), "max_imag": ev.imag.max(),
            "max_mag": np.abs(ev).max()}


def bspline_filter_operator(n: int = 64, degree: int = 7, p: int = 6, lam2=1e-6) -> np.ndarray:
    """Dense B-spline filter on [0, 2 pi] with n points (scalar lambda^2 or vector lambda)."""
    basis = Basis1D(0.0, 2.0 * np.pi, n, degree)
    lam = np.sqrt(lam2) if np.ndim(lam2) == 0 else np.asarray(lam2)
    return filter_matrix(basis.matrix(0), basis.scaled_penalty(p), lam)


def _interior(n: int, alpha: float):
    a = (93.0 + 70.0 * alpha) / 128.0
    b = (7.0 + 18.0 * alpha) / 16.0
    c = (-7.0 + 14.0 * alpha) / 32.0
    d = (1.0 - 2.0 * alpha) / 16.0
    e = (-1.0 + 2.0 * alpha) / 128.0
    erow = np.zeros(n)
    erow[0] = a
    for k, v in enumerate((b, c, d, e), start=1):
        erow[k] = erow[-k] = v / 2.0
    irow = np.zeros(n)
    irow[0], irow[1], irow[-1] = 1.0, alpha, alpha
    E = np.array([np.roll(erow, i) for i in range(n)])
    I = np.array([np.roll(irow, i) for i in range(n)])
    for M in (E, I):
        M[0, :] = 0.0
        M[-1, :] = 0.0
        M[0, 0] = M[-1, -1] = 1.0
    return E, I


def compact_filter_operator(n: int = 64, boundary: str = "8th", alpha: float = 0.49) -> np.ndarray:
    """8th-order compact filter with the named boundary closure: '8th', '6th', '4th' or '2nd-centered'."""
    E, I = _interior(n, alpha)
    al = alpha
    if boundary == "8th":
        rows = [
            [1/256 + 127*al/128, 31/32 + al/16, 7/64 + 25*al/32, -7/32 + 7*al/16,
             7*(5 - 10*al)/128, -7/32 + 7*al/16, 7/64 - 7*al/32, -1/32 + al/16, 1/256 - al/128],
            [-1/256 + al/128, 1/32 + 15*al/16, 57/64 + 7*al/32, 7/32 + 9*al/16,
             -7*(5 - 10*al)/128, 7/32 - 7*al/16, -7/64 + 7*al/32, 1/32 - al/16, -1/256 + al/128],
            [1/256 - al/128, -1/32 + al/16, 7/64 + 25*al/32, 25/32 + 7*al/16,
             35/128 + 29*al/64, -7/32 + 7*al/16, 7/64 - 7*al/32, -1/32 + al/16, 1/256 - al/128],
        ]
        first = 0
    elif boundary == "6th":
        r3 = [-1/64 + al/32, 3/32 + 13*al/16, 49/64 + 15*al/32, 5/16 + 3*al/8,
              -15/64 + 15*al/32, 3/32 - 3*al/16, -1/64 + al/32]
        rows = [[1/64 + 31*al/32, 29/32 + 3*al/16, 15/64 + 17*al/32, -5/16 + 5*al/32,
                 15/64 - 15*al/32, -3/32 + 3*al/16, 1/64 - al/32], r3, r3]
        first = 0
    elif boundary == "4th":
        r = [1/16 + 7*al/8, 3/4 + al/2, 3/8 + al/4, -1/4 + al/2, 1/16 - al/8]
        rows = [r, r, r]
        first = 0
    elif boundary == "2nd-centered":
        rows = [
            [(1/2 + al)/2, (1/2 + al), (1/2 + al)/2],
            [(-1/8 + al/4)/2, (1/2 + al)/2, (5/8 + 3*al/4), (1/2 + al)/2, (-1/8 + al/4)/2],
            [(1/32 - al/16)/2, (-3/16 + 3*al/8)/2, (15/32 + 17*al/16)/2, (11/16 + 5*al/8),
             (15/32 + 17*al/16)/2, (-3/16 + 3*al/8)/2, (1/32 - al/16)/2],
        ]
        first = 1
    else:
        raise ValueError(boundary)
    up = np.zeros((3, n))
    for k, r in enumerate(rows):
        up[k, : len(r)] = r
    E[first:first + 3, :] = up
    E[n - first - 3:n - first, :] = np.flipud(np.fliplr(up))
    return np.linalg.solve(I, E)


def stability_table(n: int = 64) -> list:
    """Spectral summaries for the B-spline filter and the compact-filter closures."""
    rows = [("8th-order B-spline (p = 6, lambda^2 = 1e-6)",
             spectrum_summary(bspline_filter_operator(n)))]
    for bnd in ("8th", "6th", "4th", "2nd-centered"):
        rows.append((f"Compact, {bnd} boundary", spectrum_summary(compact_filter_operator(n, bnd))))
    return rows
