"""Boundary-consistent B-spline filtering with GCV-selected filter strength.

Reference implementation of the filter in Bay et al., Journal of Computational
Physics 419 (2020) 109680. Modules:

  basis      clamped B-spline collocation bases and derivative penalty matrices
  filter     the filter A = B (B^T B + Omega)^{-1} B^T with LAPACK banded LU,
             Kronecker line operators in several dimensions
  gcv        generalized cross validation: Newton's method for scalar strengths,
             analytic gradient for a vector strength
  inverse    warm-started, cubically convergent iterative inverse
  stability  eigenvalue analysis on bounded grids, compact-filter closures
  burgers2d  under-resolved 2-D Burgers test of the paper
  verify     one-command reproduction of the headline numbers
"""
from .basis import Basis1D
from .filter import BandedLU, Filter1D, filter_matrix, filter_nd
from .gcv import ScalarGCV, SpectralDirection, VectorGCV
from .inverse import cubic_inverse

__all__ = ["Basis1D", "BandedLU", "Filter1D", "filter_matrix", "filter_nd",
           "ScalarGCV", "SpectralDirection", "VectorGCV", "cubic_inverse"]
