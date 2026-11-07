"""Two-dimensional viscous Burgers' equation on an under-resolved B-spline grid.

    u_t + u u_x + v u_y = nu (u_xx + u_yy),   v_t + u v_x + v v_y = nu (v_xx + v_yy)

on [0, 2 pi]^2, nu = 1/80, t in [0, 1], dt = 1e-2, with the exact traveling-front
solution supplying the initial and time-dependent Dirichlet data. Degree-7
B-spline collocation on 20 x 20 points (paper Section 5.3, dissertation Table 14.1)
is deliberately under-resolved so that aliasing errors grow; the B-spline filter
is applied once per time step with Kronecker line operators.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla

from .basis import Basis1D
from .filter import Filter1D, filter_nd

#: SSP Runge-Kutta weights (6 stages) used by the archived driver.
RK6 = np.array([1.0 / 9.0, 2.0 / 5.0, 0.0, 4.0 / 9.0, 0.0, 2.0 / 45.0])


def exact(component: str, nu: float, x, y, t):
    """Exact traveling-front solution."""
    e = np.exp((-t - 4.0 * x + 4.0 * y) / (32.0 * nu))
    return 0.75 - 0.25 / (1.0 + e) if component == "u" else 0.75 + 0.25 / (1.0 + e)


def strength_schedule(case: str, step: int) -> float:
    """Filter strength lambda^2 per step for the Table 14.1 cases.

    ``optimized-scalar`` and ``optimized-vector`` are the strengths selected for
    the two optimized cases of Table 14.1 (both directions share the value).
    """
    if case == "unfiltered" or step < 20:
        return 0.0
    if case == "optimized-vector":
        return 1.5e-5
    if case == "optimized-scalar":
        return 1.5e-5 if step <= 60 else 1.0e-4
    raise ValueError(case)


def solve(case: str = "unfiltered", degree: int = 7, n: int = 20, nu: float = 1.0 / 80.0,
          T: float = 1.0, dt: float = 1.0e-2, verbose: bool = False):
    """Run one case: 'unfiltered', 'optimized-scalar' or 'optimized-vector'.

    Returns a dict with the final L-infinity errors of u and v and the strengths used.
    """
    basis = Basis1D(0.0, 2.0 * np.pi, n, degree)
    x1 = basis.x
    B, D, L = basis.matrix(0), basis.matrix(1), basis.matrix(2)
    P = basis.penalty(2) * (2.0 * np.pi / n) ** 2

    # Nodal ordering: field[i_x, i_y] flattened in C order, so kron(Mx, My) acts on it.
    X, Y = np.meshgrid(x1, x1, indexing="ij")
    Bf = np.kron(B, B)
    Bfi = np.linalg.inv(Bf)
    Dx, Dy = np.kron(D, B), np.kron(B, D)
    Lap = np.kron(L, B) + np.kron(B, L)

    ring = np.zeros((n, n), dtype=bool)
    ring[0, :] = ring[-1, :] = ring[:, 0] = ring[:, -1] = True
    bnd = np.flatnonzero(ring.ravel())
    inn = np.flatnonzero(~ring.ravel())

    conv_x, conv_y, diff = Dx @ Bfi, Dy @ Bfi, nu * (Lap @ Bfi)
    Cx, Cy, Cd = conv_x[np.ix_(inn, bnd)], conv_y[np.ix_(inn, bnd)], diff[np.ix_(inn, bnd)]
    Bi = Bf[np.ix_(inn, inn)]
    Dxi, Dyi, Li = Dx[np.ix_(inn, inn)], Dy[np.ix_(inn, inn)], Lap[np.ix_(inn, inn)]
    lu = sla.lu_factor(Bi)

    comps = ("u", "v")
    stage = len(RK6)
    m_in = inn.size
    U = np.zeros((2, stage, m_in))
    t = 0.0
    for d, c in enumerate(comps):
        U[d, 0] = exact(c, nu, X, Y, t).ravel()[inn]

    strengths = []
    N = int(round(T / dt))
    bvals = np.zeros((2, bnd.size))
    for step in range(N):
        for m in range(stage):
            t += RK6[m] * dt
            for d, c in enumerate(comps):
                bvals[d] = exact(c, nu, X, Y, t).ravel()[bnd]
                alpha = sla.lu_solve(lu, U[d, 0])
                bc = (-U[0, m] * (Cx @ bvals[d]) - U[1, m] * (Cy @ bvals[d]) + Cd @ bvals[d])
                rhs = (U[d, m] - 0.5 * dt * U[0, m] * (Dxi @ alpha)
                       - 0.5 * dt * U[1, m] * (Dyi @ alpha) + 0.5 * dt * nu * (Li @ alpha)
                       + 0.5 * dt * bc)
                if m + 1 != stage:
                    U[d, m + 1] = rhs
                else:
                    U[d, 0] = RK6[: m - 1] @ U[d, : m - 1] + RK6[-1] * rhs

        if case != "unfiltered":
            full = np.zeros((2, n * n))
            for d in range(2):
                full[d, inn], full[d, bnd] = U[d, 0], bvals[d]
            mu = np.full(2, strength_schedule(case, step))
            strengths.append(mu.copy())
            if np.any(mu > 0.0):
                fx, fy = Filter1D(B, P, np.sqrt(mu[0])), Filter1D(B, P, np.sqrt(mu[1]))
                for d in range(2):
                    U[d, 0] = filter_nd(full[d].reshape(n, n), [fx, fy]).ravel()[inn]

        if verbose:
            print(f"t = {t:.2f}")

    linf = []
    for d, c in enumerate(comps):
        full = np.zeros(n * n)
        full[inn], full[bnd] = U[d, 0], bvals[d]
        linf.append(np.max(np.abs(exact(c, nu, X, Y, t).ravel() - full)))
    return {"case": case, "linf_u": linf[0], "linf_v": linf[1], "strengths": np.array(strengths)}
