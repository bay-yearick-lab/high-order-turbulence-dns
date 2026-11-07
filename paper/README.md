# Paper figures

Figures of Bay et al., Journal of Computational Physics 419 (2020), rendered from the original plots.

| File | Content |
| --- | --- |
| `transfer-lambda.png` | Transfer function of the 8th-order B-spline filter (6th-derivative penalty) for λ² = 1e-1 to 1e-4 |
| `transfer-deriv.png` | Transfer function at λ² = 1e-6 for 2nd-, 4th- and 6th-derivative penalties |
| `transfer-hyper-bspline.png` | B-spline filter versus 6th-derivative hyperviscosity |
| `transfer-compact-8-area-bspline.png` | B-spline filter versus 8th-order explicit and compact filters |
| `transfer-optimized-4-area-bspline.png` | B-spline filter versus an optimized 9-point finite-difference filter |
| `burgers-7.png`, `zoom-burgers-7.png` | 1-D Burgers at t = 4 with degree-7 B-splines, 32 points: no filter, constant, optimized scalar and optimized vector strength |
| `energy-removal.png` | Energy removed per step by each filter, 1-D Burgers |
| `change-burgers-7.png` | Pre- minus post-filter difference at t = 4: the optimized filters act only at the front |
| `burgers_turbulence.png` | Decaying Burgers turbulence spectra at t = 1 for 256, 1024 and 4096 points against DNS on 32768 points |
