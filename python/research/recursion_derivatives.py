import numpy as np


# Algorithm 2.5
#   Standard B-spline algorithm
#
# Compute derivatives of basis function N_{i,p}.
#
# Input: p is the degree of polynomial.
#        m is the highest index of U (m + 1) knots.
#        U is the knot vector.
#        i is the interval u's in from findspan subroutine.
#        u is a point in the set U (not necessary a knot).
#        n is the number of splines.
#
# Output: ders is array containing the basis functions.

def recursion_derivatives(p, m, U, i, u, n, ders):

    N = np.zeros((p + 1, p + 1))
    ND = np.zeros(p + 1)

    if ((i == 0 and u == U[0]) or \
        (i == m-p-1 and u == U[m])):
        for k in range(n + 1):
            ders[k] = 1.0
        return

    # local property
    if (u < U[i] or u >= U[i + p + 1]):
        for k in range(n + 1):
            ders[k] = 0.0
        return

    # initialize zeroth-degree functions
    for j in range(p + 1):
        if (u >= U[i + j] and u < U[i + j + 1]):
            N[j, 0] = 1.0
        else:
            N[j, 0] = 0.0
    # compute full triangular table
    for k in range(1, p + 1):
        if (N[0, k - 1] == 0.0):
            saved = 0.0
        else:
            saved = ((u - U[i]) * N[0, k - 1]) / (U[i + k] - U[i])
        for j in range(p - k + 1):
            Uleft = U[i + j + 1]
            Uright = U[i + j + k + 1]
            if (N[j + 1, k - 1] == 0.0):
                N[j, k] = saved
                saved = 0.0
            else:
                temp = N[j + 1, k - 1] / (Uright - Uleft)
                N[j, k] = saved + (Uright - u) * temp
                saved = (u - Uleft) * temp
    # the function value
    ders[0] = N[0, p]
#    # compute the derivatives
#    for k in range(1, n + 1):
#        # load appropriate column
#        for j in range(k + 1):
#            ND[j] = N[j, p - k]
#        # load table of width k
#        for jj in range(1, k + 1):
#            if (ND[0] == 0.0):
#                saved = 0.0
#            else:
#                saved = ND[0] / (U[i + p - k + jj] - U[i])
#            for j in range(k - jj + 1):
#                if (ND[j + 1] == 0.0):
#                    ND[j] = (p - k + jj) * saved
#                    saved = 0.0
#                else:
#                    Uleft = U[i + j + 1]
#                    Uright = U[i + j + p - k + jj + 1]
#                    temp = ND[j + 1] / (Uright - Uleft)
#                    ND[j] = (p - k + jj) * (saved - temp)
#                    saved = temp
#        # kth derivative
#        ders[k] = ND[0]
