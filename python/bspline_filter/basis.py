"""Clamped B-spline bases, collocation points, and derivative penalty matrices.

A field on [a, b] is represented as u(x) = sum_j B_j(x) alpha_j with N clamped
B-splines of degree d on a uniform breakpoint set. Collocation at N points x_i
gives the square, banded basis matrix B[i, j] = B_j(x_i).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import brentq


def clamped_knots(a: float, b: float, n_points: int, degree: int) -> np.ndarray:
    """Open (clamped) knot vector giving exactly ``n_points`` B-splines on [a, b]."""
    breaks = np.linspace(a, b, n_points - degree + 1)
    return np.concatenate([np.full(degree, a), breaks, np.full(degree, b)])


def _basis_function(knots: np.ndarray, degree: int, j: int) -> BSpline:
    n = len(knots) - degree - 1
    c = np.zeros(n)
    c[j] = 1.0
    return BSpline(knots, c, degree, extrapolate=False)


def collocation_points(knots: np.ndarray, degree: int) -> np.ndarray:
    """Collocation points for a clamped basis.

    End points are the domain ends. The ``degree - 1`` points next to each end sit
    at the maxima of the boundary-modified B-splines; interior points are the
    breakpoints (odd degree) or breakpoint midpoints (even degree).
    """
    n = len(knots) - degree - 1
    a, b = knots[0], knots[-1]
    x = np.empty(n)
    x[0], x[-1] = a, b
    near_wall = list(range(1, degree)) + list(range(n - degree, n - 1))
    for i in near_wall:
        f = _basis_function(knots, degree, i)
        df = f.derivative()
        s = np.linspace(knots[i], knots[i + degree + 1], 2001)
        k = int(np.nanargmax(np.nan_to_num(f(s))))
        k = min(max(k, 1), len(s) - 2)
        x[i] = brentq(lambda t: df(t), s[k - 1], s[k + 1], xtol=1e-15, rtol=1e-15)
    unique = np.unique(knots)
    for i in range(degree, n - degree):
        above = unique[unique > x[i - 1] + 1e-14]
        x[i] = above[0] if degree % 2 == 1 else 0.5 * (above[0] + above[1])
    return x


def basis_matrix(knots: np.ndarray, degree: int, x: np.ndarray, deriv: int = 0) -> np.ndarray:
    """Dense matrix of ``deriv``-th derivatives, M[i, j] = B_j^(deriv)(x_i)."""
    n = len(knots) - degree - 1
    out = np.zeros((len(x), n))
    for j in range(n):
        f = _basis_function(knots, degree, j)
        g = f.derivative(deriv) if deriv else f
        vals = g(x)
        out[:, j] = np.nan_to_num(vals)
    # Right end point: evaluate as a left limit so the last spline is included.
    if np.isclose(x[-1], knots[-1]):
        xe = np.array([knots[-1]])
        for j in range(n):
            f = BSpline(knots, np.eye(n)[j], degree, extrapolate=True)
            g = f.derivative(deriv) if deriv else f
            out[-1, j] = g(xe)[0]
    return out


def penalty_matrix(knots: np.ndarray, degree: int, p: int) -> np.ndarray:
    """Penalty matrix P[i, j] = int B_i^(p)(x) B_j^(p)(x) dx (exact Gauss quadrature)."""
    n = len(knots) - degree - 1
    spans = np.unique(knots)
    gx, gw = np.polynomial.legendre.leggauss(degree + 1)
    xs, ws = [], []
    for lo, hi in zip(spans[:-1], spans[1:]):
        xs.append(0.5 * (hi - lo) * gx + 0.5 * (hi + lo))
        ws.append(0.5 * (hi - lo) * gw)
    xq, wq = np.concatenate(xs), np.concatenate(ws)
    D = np.zeros((len(xq), n))
    for j in range(n):
        D[:, j] = np.nan_to_num(_basis_function(knots, degree, j).derivative(p)(xq))
    return D.T @ (wq[:, None] * D)


@dataclass
class Basis1D:
    """Clamped B-spline collocation basis on one coordinate direction."""

    a: float
    b: float
    n: int
    degree: int

    def __post_init__(self) -> None:
        self.knots = clamped_knots(self.a, self.b, self.n, self.degree)
        self.x = collocation_points(self.knots, self.degree)

    def matrix(self, deriv: int = 0) -> np.ndarray:
        """Collocation matrix of the ``deriv``-th derivative."""
        return basis_matrix(self.knots, self.degree, self.x, deriv)

    def penalty(self, p: int) -> np.ndarray:
        """Penalty matrix of the ``p``-th derivative."""
        return penalty_matrix(self.knots, self.degree, p)

    def scaled_penalty(self, p: int) -> np.ndarray:
        """Penalty matrix scaled by h^(2p-1), h the breakpoint spacing.

        With this scaling a given lambda^2 targets the same normalized wavenumber
        k h on every grid, so the filter strength is resolution independent.
        """
        h = (self.b - self.a) / (self.n - self.degree)
        return self.penalty(p) * h ** (2 * p - 1)

    @property
    def bandwidth(self) -> int:
        """Half bandwidth of the collocation matrix (``degree``)."""
        return self.degree
