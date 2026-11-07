"""Cubic iterative inverse and GCV derivatives."""
import numpy as np

from bspline_filter.basis import Basis1D
from bspline_filter.filter import filter_matrix
from bspline_filter.gcv import ScalarGCV, SpectralDirection, VectorGCV
from bspline_filter.inverse import cubic_inverse


def test_cubic_convergence_rate():
    rng = np.random.default_rng(0)
    A = rng.standard_normal((30, 30)) + 30 * np.eye(30)
    V0 = np.linalg.inv(A) @ (np.eye(30) + 1e-2 * rng.standard_normal((30, 30)) / 30)
    _, hist = cubic_inverse(A, V0)
    assert hist[-1] < 1e-13
    # E_{k+1} = E_k^3: the log residual roughly triples per iteration.
    assert np.log(hist[1]) / np.log(hist[0]) > 2.5


def test_spectral_form_equals_dense_filter():
    b = Basis1D(0.0, 2.0 * np.pi, 20, 7)
    B, P = b.matrix(0), b.scaled_penalty(6)
    d = SpectralDirection.from_operators(B, P)
    mu = 3e-4
    A = d.S @ np.diag(1.0 / (1.0 + mu * d.s)) @ d.S.T
    assert np.allclose(A, filter_matrix(B, P, np.sqrt(mu)), atol=1e-10)


def test_scalar_gcv_derivatives_and_newton():
    b = Basis1D(0.0, 2.0 * np.pi, 20, 7)
    B, P = b.matrix(0), b.scaled_penalty(6)
    X, Y = np.meshgrid(b.x, b.x, indexing="ij")
    u = np.sin(X) * np.cos(Y) + 1e-3 * np.random.default_rng(4).standard_normal(X.shape)
    g = ScalarGCV([SpectralDirection.from_operators(B, P)] * 2)
    ut2 = g.transform(u) ** 2
    th = np.log([2e-3, 5e-4])
    V, gr, H = g.value_grad_hess(th, ut2)
    e = 1e-5
    for a in range(2):
        dp = g.value_grad_hess(th + e * np.eye(2)[a], ut2)
        dm = g.value_grad_hess(th - e * np.eye(2)[a], ut2)
        assert np.isclose(gr[a], (dp[0] - dm[0]) / (2 * e), rtol=1e-5)
        assert np.allclose(H[a], (dp[1] - dm[1]) / (2 * e), rtol=1e-4, atol=1e-12)
    mu, it = g.optimize(u)
    _, gr_opt, H_opt = g.value_grad_hess(np.log(mu), ut2)
    assert it < 30
    assert np.all(np.linalg.eigvalsh(H_opt) > 0)


def test_vector_gcv_gradient():
    b = Basis1D(0.0, 2.0 * np.pi, 16, 7)
    B, P = b.matrix(0), b.scaled_penalty(6)
    u = np.exp(np.sin(b.x)) + 1e-2 * np.random.default_rng(5).standard_normal(16)
    lam = 0.05 + 0.02 * np.random.default_rng(6).random(16)
    vg = VectorGCV(B, P, lam)
    V, grad = vg.value_and_grad(lam, u)
    e = 1e-7
    for k in (0, 5, 15):
        dl = np.zeros(16)
        dl[k] = e
        fd = (vg.value_and_grad(lam + dl, u)[0] - vg.value_and_grad(lam - dl, u)[0]) / (2 * e)
        assert np.isclose(grad[k], fd, rtol=1e-4, atol=1e-12)
