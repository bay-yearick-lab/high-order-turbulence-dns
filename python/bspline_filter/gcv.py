"""Automatic filter-strength selection by generalized cross validation (GCV).

    V(lambda) = (1/N ||u - A u||^2) / (1/N tr(I - A))^2

balances the change made by the filter (numerator) against the number of
degrees of freedom it removes (denominator). Minimizing V every time step makes
the filter strength a data-driven, per-step solver control.

Scalar strength, one value mu_a = lambda_a^2 per direction
    Pre-compute, once per direction, the eigen-decomposition
    (B^T B)^{-1/2} P (B^T B)^{-1/2} = U diag(s) U^T and the orthogonal matrix
    S = B (B^T B)^{-1/2} U. Then A_a = S diag(1 / (1 + mu_a s)) S^T, so V, its
    gradient and its Hessian are closed-form sums over the transformed data
    u~ = (S_1^T x S_2^T x ...) u. No matrix is rebuilt or inverted at run time,
    and Newton's method on log(mu) converges in a few iterations.

Vector strength, one value per basis function (1-D)
    A(lambda) = B (B^T B + (lambda lambda^T) o P)^{-1} B^T. The inverse is
    refreshed with the warm-started cubic iteration of ``inverse.py``; the
    analytic gradient follows from dA/dlambda_k = -B Q^{-1} (dQ/dlambda_k) Q^{-1} B^T.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .inverse import cubic_inverse


def _apply_axis(M: np.ndarray, u: np.ndarray, axis: int) -> np.ndarray:
    return np.moveaxis(np.tensordot(M, u, axes=([1], [axis])), 0, axis)


@dataclass
class SpectralDirection:
    """Pre-computed GCV spectrum of one direction: A(mu) = S diag(1/(1+mu s)) S^T."""

    S: np.ndarray
    s: np.ndarray

    @classmethod
    def from_operators(cls, B: np.ndarray, P: np.ndarray) -> "SpectralDirection":
        """Build the spectrum from the collocation matrix B and penalty matrix P."""
        w, Z = np.linalg.eigh(B.T @ B)
        Mih = (Z / np.sqrt(w)) @ Z.T
        s, U = np.linalg.eigh(Mih @ P @ Mih)
        s = np.clip(s, 0.0, None)
        return cls(S=B @ Mih @ U, s=s)


@dataclass
class ScalarGCV:
    """GCV objective with one scalar filter strength per direction."""

    dirs: list
    log_mu_bounds: tuple = (np.log(1e-14), np.log(1e2))
    history: list = field(default_factory=list)

    def transform(self, u: np.ndarray) -> np.ndarray:
        """u~ = (S_1^T x S_2^T x ...) u."""
        ut = np.asarray(u, dtype=float)
        for a, d in enumerate(self.dirs):
            ut = _apply_axis(d.S.T, ut, a)
        return ut

    def _factors(self, theta):
        shape = [len(d.s) for d in self.dirs]
        sigma = np.ones(shape)
        g, dg = [], []
        for a, (d, th) in enumerate(zip(self.dirs, theta)):
            mu = np.exp(th)
            q = mu * d.s / (1.0 + mu * d.s)
            view = [1] * len(shape)
            view[a] = -1
            sigma = sigma * (1.0 / (1.0 + mu * d.s)).reshape(view)
            g.append((-q).reshape(view))
            dg.append((-q * (1.0 - q)).reshape(view))
        return sigma, g, dg

    def value_grad_hess(self, theta, ut2):
        """V, dV/dtheta and d2V/dtheta2 at theta = log(mu), with ut2 = u~^2."""
        m = len(self.dirs)
        N = ut2.size
        sigma, g, dg = self._factors(theta)
        one = 1.0 - sigma
        d1 = [sigma * g[a] for a in range(m)]
        V1 = np.sum(one**2 * ut2) / N
        V2 = (N - np.sum(sigma)) / N
        V1a = np.array([np.sum(-2.0 * one * d1[a] * ut2) / N for a in range(m)])
        V2a = np.array([-np.sum(d1[a]) / N for a in range(m)])
        V1ab = np.zeros((m, m))
        V2ab = np.zeros((m, m))
        for a in range(m):
            for b in range(m):
                d2 = sigma * g[a] * g[b] + (sigma * dg[a] if a == b else 0.0)
                V1ab[a, b] = np.sum((2.0 * d1[a] * d1[b] - 2.0 * one * d2) * ut2) / N
                V2ab[a, b] = -np.sum(d2) / N
        V = V1 / V2**2
        grad = V1a / V2**2 - 2.0 * V1 * V2a / V2**3
        hess = (V1ab / V2**2
                - 2.0 * (np.outer(V1a, V2a) + np.outer(V2a, V1a)) / V2**3
                - 2.0 * V1 * V2ab / V2**3
                + 6.0 * V1 * np.outer(V2a, V2a) / V2**4)
        return V, grad, hess

    def optimize(self, u: np.ndarray, mu0=None, gtol: float = 1e-8, max_iter: int = 50):
        """Newton's method on theta = log(mu) with an eigenvalue-shifted Hessian
        and Armijo backtracking. Stops when ||dV/dtheta||_inf <= gtol * V.

        Returns (mu, n_iterations).
        """
        ut2 = self.transform(u) ** 2
        lo, hi = self.log_mu_bounds
        theta = np.log(np.full(len(self.dirs), 1e-6) if mu0 is None else np.asarray(mu0, float))
        theta = np.clip(theta, lo, hi)
        self.history = []
        it = 0
        for it in range(1, max_iter + 1):
            V, gr, H = self.value_grad_hess(theta, ut2)
            self.history.append(V)
            free = ~(((theta <= lo) & (gr > 0)) | ((theta >= hi) & (gr < 0)))
            if np.max(np.abs(gr[free]), initial=0.0) <= gtol * V:
                break
            w, Z = np.linalg.eigh(H)
            w = np.maximum(np.abs(w), 1e-12 * max(np.abs(w).max(), 1e-300))
            step = -(Z / w) @ (Z.T @ gr)
            step[~free] = 0.0
            t = 1.0
            while t > 1e-12:
                trial = np.clip(theta + t * step, lo, hi)
                if self.value_grad_hess(trial, ut2)[0] <= V + 1e-4 * t * (gr @ step):
                    break
                t *= 0.5
            new = np.clip(theta + t * step, lo, hi)
            if np.max(np.abs(new - theta)) < 1e-12:
                break
            theta = new
        return np.exp(theta), it


class VectorGCV:
    """GCV objective for a vector filter strength in one direction.

    Keeps Q^{-1} between calls and refreshes it with the warm-started cubic
    iteration; ``inverse_iterations`` records the iterations used per refresh.
    """

    def __init__(self, B: np.ndarray, P: np.ndarray, lam0: np.ndarray, tol: float = 1e-14):
        self.B, self.P, self.tol = B, P, tol
        self.BtB = B.T @ B
        self.lam = np.asarray(lam0, dtype=float).copy()
        self.Qinv = np.linalg.inv(self._Q(self.lam))
        self.inverse_iterations = []

    def _Q(self, lam):
        return self.BtB + np.outer(lam, lam) * self.P

    def refresh(self, lam):
        """Update Q^{-1} for a new ``lam`` starting from the previous inverse."""
        lam = np.asarray(lam, dtype=float)
        self.Qinv, hist = cubic_inverse(self._Q(lam), self.Qinv, tol=self.tol)
        self.inverse_iterations.append(len(hist) - 1)
        self.lam = lam.copy()

    def value_and_grad(self, lam, u):
        """V(lam) and its analytic gradient."""
        self.refresh(lam)
        N = u.size
        B, P, Qi = self.B, self.P, self.Qinv
        c = Qi @ (B.T @ u)
        r = u - B @ c
        a = Qi @ (B.T @ r)
        V1 = r @ r / N
        trA = np.sum(Qi * self.BtB)
        V2 = (N - trA) / N
        dV1 = (2.0 / N) * (a * (P @ (lam * c)) + c * (P @ (lam * a)))
        G = Qi @ self.BtB @ Qi
        dV2 = (2.0 / N) * ((G * P) @ lam)
        V = V1 / V2**2
        return V, dV1 / V2**2 - 2.0 * V1 * dV2 / V2**3

    def optimize(self, u, lam_max: float = 1e-2, max_iter: int = 200):
        """Minimize V over 0 <= lam_k <= lam_max with the analytic gradient (L-BFGS-B)."""
        res = minimize(self.value_and_grad, self.lam, args=(u,), jac=True, method="L-BFGS-B",
                       bounds=[(0.0, lam_max)] * len(self.lam), options={"maxiter": max_iter})
        self.refresh(res.x)
        return res.x, res
