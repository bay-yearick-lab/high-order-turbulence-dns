"""Basis and penalty-matrix identities."""
import numpy as np
import pytest

from bspline_filter.basis import Basis1D


@pytest.mark.parametrize("degree", [3, 5, 7])
def test_partition_of_unity_and_band(degree):
    b = Basis1D(0.0, 2.0 * np.pi, 32, degree)
    B = b.matrix(0)
    assert np.allclose(B.sum(axis=1), 1.0, atol=1e-14)
    i, j = np.nonzero(np.abs(B) > 0)
    assert np.max(np.abs(i - j)) <= degree


@pytest.mark.parametrize("p", [2, 4, 6])
def test_penalty_symmetric_psd_with_polynomial_null_space(p):
    b = Basis1D(-1.0, 1.0, 16, 7)
    P = b.penalty(p)
    scale = np.abs(P).max()
    assert np.allclose(P, P.T, rtol=0, atol=1e-12 * scale)
    # Polynomials of degree < p have zero p-th derivative, so they are not penalized.
    B = b.matrix(0)
    for k in range(p + 1):
        c = np.linalg.solve(B, b.x**k)
        energy = c @ P @ c
        if k < p:
            assert abs(energy) < 1e-9 * scale * (c @ c)
        else:
            exact = 2.0 * float(np.prod(np.arange(1, p + 1))) ** 2  # int_{-1}^{1} (p!)^2 dx
            assert energy == pytest.approx(exact, rel=1e-6)


def test_derivative_matrix_exact_for_polynomials():
    b = Basis1D(0.0, 2.0, 24, 7)
    coeffs = np.linalg.solve(b.matrix(0), b.x**5)
    assert np.allclose(b.matrix(1) @ coeffs, 5 * b.x**4, atol=1e-10)
    assert np.allclose(b.matrix(2) @ coeffs, 20 * b.x**3, atol=1e-8)
