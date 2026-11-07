import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from findspan import findspan
from derivatives import derivatives
from knot_multiplicity import knot_multiplicity


def B_idx(u, N, x, n, n_splines, U, degree):

    i = findspan(n_splines - 1, degree, u, U)
    derivatives(i, u, degree, n, U, x)

    a = (np.abs(x[0, :] - N)).argmin()
    return a


def nonperiodic_init(mink, maxk, points, degree, deriv, U2=None):

    pts = points - degree + 1

    U = np.linspace(mink, maxk, pts)
    # knot multiplicity
    U = knot_multiplicity(degree, U)

    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # collocation pts
    pts = n_splines
    # collocation points
    co = np.zeros((2, pts))

    collocation_points(co, U, degree)
    # random u value whom we are trying to find the ith interval it lies in
    u = co[0, :]

    if U2 is None:
        U2 = U.copy()

    n_splines2 = len(U2) - degree - 1

    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((pts, deriv + 1, n_splines2))

    for j in range(len(u)):
        i = findspan(n_splines2 - 1, degree, u[j], U2)
        derivatives(i, u[j], degree, deriv, U2, ders[j, :, (i - degree):(degree + i)])

    return co[0, :], ders, U2


def collocation_points(co, U, degree, string=None):
    from galerkin import Bi
    import scipy.optimize as so

    if (string == 'periodic'):
        pts = len(co[0, :])
        for i in range(pts):
            if (i == 0):
                co[0, i] = U[0]
            elif (i == pts - 1):
                co[0, i] = U[-1]
            elif (degree%2 == 0):
                co[0, i] = 0.5 * (U[i-1] + U[i])
            else:
                co[0, i] = U[i]

    else:
        pts = len(co[0, :])
        xx = 20 * pts
        # number of B-spline curves
        n_splines = len(U) - degree - 1
        u = np.linspace(min(U), max(U), xx)
        # der[i, j] is the kth derivative where
        # i is the number of distinct knots and j is the degree
        N = np.zeros((xx, 1, n_splines))
        # temporary array
        x = np.zeros((2, degree + 1))

        for j in range(xx):
            i = findspan(n_splines - 1, degree, u[j], U)
            derivatives(
                i, u[j], degree, 0, U, N[j, :, (i - degree):(degree + i)])

        pts = len(co[0, :])
        for i in range(pts):
            if (i == 0):
                co[0, i] = u[0]
            elif ((1 <= i < degree) or ((pts - degree) <= i < (pts - 1))):
                idx = np.argmax(N[:, 0, i])
                a = B_idx(u[idx], N[idx, 0, i], x, 1, n_splines, U, degree)
                co[0, i] = so.newton(Bi, u[idx], args=(a, x, 1, n_splines, U, degree), tol=1e-12)
            elif (i == pts - 1):
                co[0, i] = u[-1]
            else:
                value = np.where(U > co[0, i-1])[0]
                if (degree%2 == 0):
                    co[0, i] = 0.5 * (U[value[0]] + U[value[1]])
                else:
                    co[0, i] = U[value[0]]
    return


def bvp_boundary_layer(pts):
    # degree
    degree = 5
    mink = -1.0
    maxk = 1.0
    knotpts = pts
    # knot vector
    U = np.linspace(mink, maxk, knotpts)
    U = knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # collocation pts
    pts = n_splines
    # collocation points
    co = np.zeros((2, pts))
    # number of derivatives
    n = 2

    xx = 100000
    u = np.linspace(mink, maxk, xx)
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((xx, n + 1, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    collocation_points(co, u, U, ders, degree)

    u = co[0, :]
    ders = np.zeros((pts, n + 1, n_splines))
    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    # bc
    bc0 = 1.0
    bc1 = 0.0

    M = ders[:, 2, :] - 100.0 * ders[:, 0, :]

    # delete first/last row
    M = np.delete(M, 0, 0)
    M = np.delete(M, -1, 0)

    # save first/last column
    firstcolumn = M[:, 0]
    lastcolumn = M[:, -1]

    # delete first/last column
    M = np.delete(M, 0, 1)
    M = np.delete(M, -1, 1)

    n, m = M.shape
    assert(m == n)

    b = -100.0 * np.ones(n)

    # place bc to rhs
    b = b - bc0 * firstcolumn
    b = b - bc1 * lastcolumn

    alpha = la.solve(M, b)

    # put back bc
    alpha = np.insert(alpha, 0 , bc0)
    alpha = np.append(alpha, bc1)

    # interpolation error
    xx = 10000
    u = np.linspace(mink, maxk, xx)
    ders = np.zeros((len(u), n + 1, n_splines))
    y = np.zeros(len(u))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(i + degree)])

    for j in range(len(u)):
        for i in range(len(alpha)):
            y[j] = y[j] + alpha[i] * ders[j, 0, i]

    sol = np.ones(xx) - np.sinh(10.0 * (u + 1.0)) / np.sinh(2.0 * 10.0)
    error = la.norm(sol - y, np.inf)
    print(error)

#    plt.figure(1, figsize=(11, 8))
#    plt.plot(u, y, 'o')
#    plt.plot(u, sol, linewidth=3)
#    plt.show()

    return error


def bvp_smooth_2(pts):
    # degree
    degree = 5
    mink = 0.0
    maxk = 1.0
    knotpts = pts
    # knot vector
    U = np.linspace(mink, maxk, knotpts)
    U = knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # collocation pts
    pts = n_splines
    # collocation points
    co = np.zeros((2, pts))
    # number of derivatives
    n = 2

    xx = 100000
    u = np.linspace(mink, maxk, xx)
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((xx, n + 1, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    collocation_points(co, u, U, ders, degree)

    print(U)
    print(co[0, :])
    u = co[0, :]
    ders = np.zeros((pts, n + 1, n_splines))
    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    # bc
    bc0 = 0.0
    bc1 = 0.0

    M = ders[:, 2, :] - 4.0 * ders[:, 0, :]

    # delete first/last row
    M = np.delete(M, 0, 0)
    M = np.delete(M, -1, 0)

    # save first/last column
    firstcolumn = M[:, 0]
    lastcolumn = M[:, -1]

    # delete first/last column
    M = np.delete(M, 0, 1)
    M = np.delete(M, -1, 1)

    n, m = M.shape
    assert(m == n)

    b = 4.0 * np.cosh(np.ones(n))

    # place bc to rhs
    b = b - bc0 * firstcolumn
    b = b - bc1 * lastcolumn

    alpha = la.solve(M, b)

    # put back bc
    alpha = np.insert(alpha, 0 , bc0)
    alpha = np.append(alpha, bc1)

    # interpolation error
    xx = 10000
    u = np.linspace(mink, maxk, xx)
    ders = np.zeros((len(u), n + 1, n_splines))
    y = np.zeros(len(u))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(i + degree)])

    for j in range(len(u)):
        for i in range(len(alpha)):
            y[j] = y[j] + alpha[i] * ders[j, 0, i]

    sol = np.cosh(2.0 * u - 1.0) - np.cosh(1.0)
    error = la.norm(sol - y, np.inf)
    print(error)

#    plt.figure(1, figsize=(11, 8))
#    plt.plot(u, y, 'o')
#    plt.plot(u, sol, linewidth=3)
#    plt.show()

    return error


def bvp_smooth_1(pts):
    # degree
    degree = 5
    mink = 0.0
    maxk = 1.0
    knotpts = pts
    # knot vector
    U = np.linspace(mink, maxk, knotpts)
    U = knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # collocation pts
    pts = n_splines
    # collocation points
    co = np.zeros((2, pts))
    # number of derivatives
    n = 2

    xx = 100000
    u = np.linspace(mink, maxk, xx)
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((xx, n + 1, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    collocation_points(co, u, U, ders, degree)

    u = co[0, :]
    ders = np.zeros((pts, n + 1, n_splines))
    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    # bc
    bc0 = 1.0
    bc1 = 1.0

    M = ders[:, 2, :] - 100.0 * ders[:, 0, :]

    # delete first/last row
    M = np.delete(M, 0, 0)
    M = np.delete(M, -1, 0)

    # save first/last column
    firstcolumn = M[:, 0]
    lastcolumn = M[:, -1]

    # delete first/last column
    M = np.delete(M, 0, 1)
    M = np.delete(M, -1, 1)

    n, m = M.shape
    assert(m == n)

    b = np.zeros(n)

    # place bc to rhs
    b = b - bc0 * firstcolumn
    b = b - bc1 * lastcolumn

    alpha = la.solve(M, b)

    # put back bc
    alpha = np.insert(alpha, 0 , 1.0)
    alpha = np.append(alpha, 1.0)

    # interpolation error
    xx = 10000
    u = np.linspace(mink, maxk, xx)
    ders = np.zeros((xx, n + 1, n_splines))
    y = np.zeros(xx)

    for j in range(xx):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(i + degree)])

    for j in range(xx):
        for i in range(len(alpha)):
            y[j] = y[j] + alpha[i] * ders[j, 0, i]

    sol = np.cosh(10.0 * (u - 0.5)) / np.cosh(5.0)
    error = la.norm(sol - y, np.inf)
    print(error)

#    plt.figure(1, figsize=(11, 8))
#    plt.plot(u, y, 'o')
#    plt.plot(u, sol, linewidth=3)
#    plt.show()

    return error


def test_bvp_degree_five(problem):
    import csv

    error = []
    pts = np.array([10, 20, 40, 80, 100, 150, 200, 250, 300, 350])

    for i in range(len(pts)):
        if (problem == 1):
            error.append(bvp_smooth_1(pts[i]))
        elif (problem == 2):
            error.append(bvp_smooth_2(pts[i]))
        elif (problem == 3):
            error.append(bvp_boundary_layer(pts[i]))

    temp = zip(pts, error)

    with open('output.csv', "w") as output:
        writer = csv.writer(output, lineterminator='\n', delimiter=',')
        for val in temp:
            writer.writerow(val)


def main():
    #test_bvp_degree_five(3)
    bvp_smooth_2(10)


if __name__ == "__main__":
    main()
