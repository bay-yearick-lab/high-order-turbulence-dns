"""Warm-started, cubically convergent iterative inverse.

For the vector filter strength the normal matrix Q_n = B^T B + (lam lam^T) o P
changes every time step, but only slightly. Starting from the previous inverse
V_0 = Q_{n-1}^{-1}, the third-order hyperpower iteration

    V_{k+1} = V_k (3 I - 3 Q V_k + (Q V_k)^2)

converges cubically: the residual E_k = I - Q V_k obeys E_{k+1} = E_k^3, so a
warm start with ||E_0|| = 1e-5 reaches round-off in two iterations.
"""
from __future__ import annotations

import numpy as np


def cubic_inverse(Q: np.ndarray, V0: np.ndarray, tol: float = 1e-14, max_iter: int = 20):
    """Refine an approximate inverse ``V0`` of ``Q``.

    Iterates until ||I - Q V||_inf <= tol or the residual stops decreasing
    (round-off floor). If the warm start is not contractive (residual >= 1) the
    inverse is recomputed directly. Returns ``(V, history)`` where
    ``history[k]`` is ||I - Q V_k||_inf and history[0] is the warm-start residual.
    """
    n = Q.shape[0]
    eye = np.eye(n)
    QV = Q @ V0
    history = [np.linalg.norm(eye - QV, np.inf)]
    if not history[0] < 1.0:
        V = np.linalg.inv(Q)
        return V, history + [np.linalg.norm(eye - Q @ V, np.inf)]
    V = V0.copy()
    for _ in range(max_iter):
        if history[-1] <= tol:
            break
        V_new = V @ (3.0 * eye - 3.0 * QV + QV @ QV)
        QV_new = Q @ V_new
        r = np.linalg.norm(eye - QV_new, np.inf)
        if len(history) > 1 and r > 0.5 * history[-1]:
            break
        V, QV = V_new, QV_new
        history.append(r)
    return V, history
