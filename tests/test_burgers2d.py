"""Reproduce the under-resolved 2-D Burgers errors of paper Table 9."""
import pytest

from bspline_filter import burgers2d


@pytest.mark.parametrize("case,expected", [("unfiltered", 2.53e-1),
                                           ("optimized-scalar", 2.07e-2),
                                           ("optimized-vector", 1.57e-2)])
def test_table9_errors(case, expected):
    r = burgers2d.solve(case)
    assert r["linf_u"] == pytest.approx(expected, rel=5e-3)
    assert r["linf_v"] == pytest.approx(expected, rel=5e-3)
