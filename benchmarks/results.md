# Benchmarks

Run with `make bench` (`python benchmarks/run_benchmarks.py`). Results below: Intel Core Ultra 7 265F, single process, Python 3.14, NumPy 2.5, SciPy 1.18.

| Benchmark | Time | Notes |
| --- | --- | --- |
| Coefficient solve B^-1 u, 60x80x60 (3 banded line sweeps) | 42.1 ms |  |
| B-spline filter A u, 60x80x60 (3 banded line sweeps) | 43.9 ms | 1.04x the coefficient solve |
| Newton-optimized GCV, 3 strengths, 60x80x60 | 13 Newton iterations | lambda^2 = 4.6e+00, 1.3e+01, 7.1e+01 |
| Warm-started cubic inverse, N = 80, 200 steps | 0.69 ms/step | mean 2.00, max 3 iterations |

The filter costs about as much as the coefficient solve that every collocation right-hand side already performs, because both are three sweeps of banded LAPACK solves along grid lines. In the Fortran channel solver (60 x 80 x 60, paper Table 12) the measured overhead per Runge-Kutta step is 0.23% for a constant strength, 8.7% with per-step scalar optimization and 21.9% with the wall-normal vector strength.
