import numpy as np


# Algorithm 2.3
#   Standard B-spline algorithm
#
# Computes the non-zero basis functions and their derivatives.
#
# Input: i is the interval u's in from findspan subroutine.
#        u is a point in the set U (not necessary a knot).
#        p is the degree of polynomial.
#        n is the number of splines.
#        U is the knot vector.
#        ders is array containing the basis functions.
def derivatives(i, u, p, n, U, ders):

    # create temporary arrays
    ndu = np.zeros((p + 1, p + 1))
    a = np.zeros((2, p + 1))
    left = np.zeros(p + 1)
    right = np.zeros(p + 1)

    ndu[0, 0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u

        saved = 0.0

        for r in range(j):
            # lower triangle
            ndu[j, r] = right[r + 1] + left[j - r]
            temp = ndu[r, j - 1] / ndu[j, r]

            # upper triangle
            ndu[r, j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp

        ndu[j, j] = saved

    # load basis functions
    for j in range(p + 1):
        ders[0, j] = ndu[j, p]

    # computes derivatives with the standard derivative recursion.
    #
    #   N^(k)_{i, p} (u) = p (N^(k - 1)_{i, p - 1} / (u_{i + p} - u_i)
    #                       - N^(k - 1)_{i + 1, p - 1} / (u_{i + p + 1} - u_{i + 1}))
    # loop over function index
    for r in range(p + 1):
        s1 = 0
        s2 = 1
        a[0, 0] = 1.0

        for k in range(1, n + 1):
            d = 0.0
            rk = r - k
            pk = p - k
            if (r >= k):
                a[s2, 0] = a[s1, 0] / ndu[pk + 1, rk]
                d = a[s2, 0] * ndu[rk, pk]
            if (rk >= -1):
                j1 = 1
            else:
                j1 = -rk
            if (r - 1 <= pk):
                j2 = k - 1
            else:
                j2 = p - r

            for j in range(j1, j2 + 1):
                a[s2, j] = (a[s1, j] - a[s1, j - 1]) / ndu[pk + 1, rk + j]
                d += a[s2, j] * ndu[rk + j, pk]

            if (r <= pk):
                a[s2, k] = -a[s1, k - 1] / ndu[pk + 1, r]
                d += a[s2, k] * ndu[r, pk]

            ders[k, r] = d
            j = s1
            s1 = s2
            s2 = j
    # multiply by correct factors in Eq. (2.9)
    r = p
    for k in range(1, n + 1):
        for j in range(p + 1):
            ders[k, j] *= r
        r *= (p - k)
