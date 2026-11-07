"""Boundary-consistent B-spline filter.

The filtered field minimizes ||u - B alpha||^2 + alpha^T Omega alpha, with
Omega = lambda^2 P (scalar strength) or Omega = L P L, L = diag(lambda_vec)
(vector strength, one value per basis function). The filter operator is

    A = B (B^T B + Omega)^{-1} B^T,

a symmetric positive semi-definite matrix with every eigenvalue in [0, 1]:
the filter can never amplify a mode, on any domain, for any lambda or P.
The normal matrix Q = B^T B + Omega is banded (half bandwidth = degree), so it
is factored once with LAPACK banded LU (dgbtrf) and applied with dgbtrs.
In several dimensions the filter is a Kronecker product of 1-D operators and
is applied line by line along each axis.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.linalg import lapack


def to_banded(Q: np.ndarray, kl: int, ku: int) -> np.ndarray:
    """Pack a square matrix in LAPACK general-band storage with room for LU fill."""
    n = Q.shape[0]
    ab = np.zeros((2 * kl + ku + 1, n))
    for j in range(n):
        for i in range(max(0, j - ku), min(n, j + kl + 1)):
            ab[kl + ku + i - j, j] = Q[i, j]
    return ab


def half_bandwidth(Q: np.ndarray, tol: float = 0.0) -> int:
    """Largest |i - j| with |Q[i, j]| > tol."""
    i, j = np.nonzero(np.abs(Q) > tol)
    return int(np.max(np.abs(i - j))) if i.size else 0


class BandedLU:
    """LAPACK banded LU factorization (dgbtrf) with dgbtrs solves."""

    def __init__(self, Q: np.ndarray):
        self.n = Q.shape[0]
        self.k = half_bandwidth(Q, 1e-300)
        ab = to_banded(Q, self.k, self.k)
        self.lu, self.piv, info = lapack.dgbtrf(ab, self.k, self.k)
        if info != 0:
            raise np.linalg.LinAlgError(f"dgbtrf failed, info={info}")

    def solve(self, rhs: np.ndarray, block: int = 128) -> np.ndarray:
        """Solve Q x = rhs for one or several right-hand sides (columns).

        Right-hand sides are processed in cache-sized column blocks.
        """
        rhs = np.asarray(rhs, dtype=float)
        if rhs.ndim == 1:
            rhs = rhs[:, None]
            squeeze = True
        else:
            squeeze = False
        out = np.empty_like(rhs, order="F")
        for j in range(0, rhs.shape[1], block):
            x, info = lapack.dgbtrs(self.lu, self.k, self.k,
                                    np.asfortranarray(rhs[:, j:j + block]), self.piv)
            if info != 0:
                raise np.linalg.LinAlgError(f"dgbtrs failed, info={info}")
            out[:, j:j + block] = x
        return out[:, 0] if squeeze else out


def normal_matrix(B: np.ndarray, P: np.ndarray, lam) -> np.ndarray:
    """Q = B^T B + lambda^2 P (scalar) or B^T B + diag(lam) P diag(lam) (vector)."""
    lam = np.asarray(lam, dtype=float)
    if lam.ndim == 0:
        return B.T @ B + lam**2 * P
    return B.T @ B + np.outer(lam, lam) * P


def filter_matrix(B: np.ndarray, P: np.ndarray, lam) -> np.ndarray:
    """Dense filter operator A = B Q^{-1} B^T (for analysis and testing)."""
    Q = normal_matrix(B, P, lam)
    return B @ np.linalg.solve(Q, B.T)


class Filter1D:
    """1-D B-spline filter with a pre-factored banded normal matrix."""

    def __init__(self, B: np.ndarray, P: np.ndarray, lam):
        self.B = B
        self.P = P
        self._Bs = sparse.csr_matrix(B)
        self._BTs = sparse.csr_matrix(B.T)
        self.set_strength(lam)

    def set_strength(self, lam) -> None:
        """Refactor for a new scalar or vector filter strength."""
        self.lam = np.asarray(lam, dtype=float)
        self.lu = BandedLU(normal_matrix(self.B, self.P, self.lam))

    def apply(self, u: np.ndarray, axis: int = 0) -> np.ndarray:
        """Filter nodal values ``u`` along ``axis`` (all other axes are lines)."""
        v = np.moveaxis(np.asarray(u, dtype=float), axis, 0)
        shape = v.shape
        rhs = self._BTs @ v.reshape(shape[0], -1)
        alpha = self.lu.solve(rhs)
        out = (self._Bs @ alpha).reshape(shape)
        return np.moveaxis(out, 0, axis)


def filter_nd(u: np.ndarray, filters) -> np.ndarray:
    """Apply one Filter1D per axis: (A_1 x A_2 x ...) u with Kronecker line operators."""
    out = np.asarray(u, dtype=float)
    for axis, f in enumerate(filters):
        if f is not None:
            out = f.apply(out, axis=axis)
    return out
