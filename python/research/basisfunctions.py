import numpy as np


# Algorithm 2.2
#   Standard B-spline algorithm
#
# Computes the non-vanishing basis functions N for a given p,
#   N[0] ... N[p].
#
# Input: i is the interval u's in from findspan subroutine.
#        u is a point in the set U (not necessary a knot).
#        p is the degree of polynomial.
#        U is the knot vector.
#        N is array containing the basis functions.
def basisfunctions(i, u, p, U, N):
    N[0] = 1.0

    # create arrays left and right
    left = np.zeros(p + 1)
    right = np.zeros(p + 1)

    for j in range(1, p + 1):
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u
        saved = 0.0

        for r in range(j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp

        N[j] = saved
