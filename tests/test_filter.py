"""Filter correctness: banded LAPACK path, polynomial preservation, Kronecker lines, order."""
import numpy as np

from bspline_filter.basis import Basis1D
from bspline_filter.filter import Filter1D, filter_matrix, filter_nd


def _setup(n=40, degree=7, p=6):
    b = Basis1D(0.0, 2.0 * np.pi, n, degree)
    return b, b.matrix(0), b.scaled_penalty(p)


def test_banded_lu_matches_dense():
    b, B, P = _setup()
    u = np.random.default_rng(0).standard_normal(b.n)
    for lam in (np.sqrt(1e-4), 1e-2 * np.ones(b.n)):
        assert np.allclose(Filter1D(B, P, lam).apply(u), filter_matrix(B, P, lam) @ u, atol=1e-11)


def test_polynomials_below_penalty_order_pass_unchanged():
    b, B, P = _setup(p=6)
    f = Filter1D(B, P, np.sqrt(1e-1))
    for k in range(6):
        u = (b.x / np.pi - 1.0) ** k
        assert np.allclose(f.apply(u), u, atol=1e-9)


def test_kronecker_line_operators():
    b, B, P = _setup(n=16)
    fx, fy = Filter1D(B, P, 1e-2), Filter1D(B, P, 3e-2)
    U = np.random.default_rng(1).standard_normal((16, 16))
    K = np.kron(filter_matrix(B, P, 1e-2), filter_matrix(B, P, 3e-2))
    assert np.allclose(filter_nd(U, [fx, fy]).ravel(), K @ U.ravel(), atol=1e-11)


def test_filter_preserves_high_order_consistency():
    ns = [24, 32, 48, 64]
    errs = []
    for n in ns:
        b, B, P = _setup(n=n)
        u = np.exp(np.sin(b.x))
        errs.append(np.max(np.abs(Filter1D(B, P, 1e-3).apply(u) - u)))
    rates = -np.diff(np.log(errs)) / np.diff(np.log(ns))
    assert rates.min() > 7.0
