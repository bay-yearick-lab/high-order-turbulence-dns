import numpy as np


# Algorithm 2.1
#   Standard B-spline algorithm
#
# Determines the knot span index, i.
#   Input: n is the number of curves, where
#               n = m - p - 1 and m is number of knots.
#          p is the degree of polynomial.
#          u is a distinct knot which is in the set U,
#               where u \in (u_{m - p + 1}, u_{m - p}].
#          U is the knot vector.
#
#
def findspan(n, p, u, U):

    # special case
    if (u == U[n + 1]):
        return(n)

    # binary search
    low = p
    high = n + 1
    mid = int((low + high) / 2)

    while (u < U[mid] or u >= U[mid + 1]):
        if (u < U[mid]):
            high = mid
        else:
            low = mid
        mid = int((low + high) / 2)
    return mid
