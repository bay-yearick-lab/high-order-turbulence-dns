"""Linear stability on a bounded grid."""
import numpy as np
import pytest

from bspline_filter.basis import Basis1D
from bspline_filter.filter import filter_matrix
from bspline_filter.stability import (bspline_filter_eigenvalues, compact_filter_operator,
                                      penalty_pencil_eigenvalues, spectrum_summary)


@pytest.mark.parametrize("n", [32, 64])
@pytest.mark.parametrize("lam2", [1e-8, 1e-4, 1e-1])
def test_bspline_filter_spectrum_in_unit_interval(n, lam2):
    b = Basis1D(0.0, 2.0 * np.pi, n, 7)
    B, P = b.matrix(0), b.scaled_penalty(6)
    nu = penalty_pencil_eigenvalues(B, P, np.sqrt(lam2))
    assert nu.min() >= -1e-14 * nu.max()
    ev = bspline_filter_eigenvalues(B, P, np.sqrt(lam2))
    assert ev.min() > 0.0 and ev.max() - 1.0 <= 1e-14 * nu.max()
    explicit = np.linalg.eigvals(filter_matrix(B, P, np.sqrt(lam2)))
    assert np.abs(explicit.imag).max() < 1e-10


def test_vector_strength_spectrum_in_unit_interval():
    b = Basis1D(0.0, 2.0 * np.pi, 48, 7)
    B, P = b.matrix(0), b.scaled_penalty(6)
    lam = 0.3 * np.random.default_rng(3).random(48)
    nu = penalty_pencil_eigenvalues(B, P, lam)
    assert nu.min() >= -1e-14 * nu.max()


@pytest.mark.parametrize("closure,expected", [("8th", 1.000426), ("6th", 1.000213)])
def test_compact_closures_amplify(closure, expected):
    s = spectrum_summary(compact_filter_operator(64, closure))
    assert s["max_mag"] == pytest.approx(expected, abs=1e-6)
    assert s["max_imag"] > 0.05


def test_low_order_centered_closure_is_contractive():
    s = spectrum_summary(compact_filter_operator(64, "2nd-centered"))
    assert s["max_mag"] <= 1.0 + 1e-12
