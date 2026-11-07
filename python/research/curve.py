import numpy as np


# Algorithm 3.1
#   Standard B-spline algorithm
#
# Computes a point on B-spline curve at a fixed u value.
#   1. find the knot span in which u lies (Algorithm A2.1)
#   2. compute the nonzero basis functions (Algorithm A2.3)
#   3. multiply the values of the nonzero basis functions with
#      corresponding control points
#
# Input:
#        n is the number of splines.
#        p is the degree of polynomial.
#        u is a point in the set U (not necessary a knot).
#        U is the knot vector.
#        P is the control point vector.
#        N is array containing the basis functions.
#
# Output:
#        C is the value of B-spline curve at point u.
#
def curvepoint(n, p, U, P, u, C, N):
    from findspan import findspan
    from basisfunctions import basisfunctions

    i = findspan(n - 1, p, u, U)
    basisfunctions(i, u, p, U, N)

    C.fill(0)
    for j in range(p + 1):
        C += N[j] * P[:, i - p + j]


def nurbs(n, p, U, P, u, C, N):
    from findspan import findspan
    from basisfunctions import basisfunctions

    i = findspan(n - 1, p, u, U)
    basisfunctions(i, u, p, U, N)

    C.fill(0)
    for j in range(p + 1):
        C += N[j] * P[:, i - p + j]
        N[j] = N[j] * P[-1, i - p + j]

    for j in range(p + 1):
        N[j] = N[j] / C[-1]

    C[0] = C[0] / C[-1]
    C[1] = C[1] / C[-1]
