import numpy as np
import numpy.linalg as la
import scipy.linalg as sl
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from findspan import findspan
from scipy.integrate import quad
from derivatives import derivatives
from basisfunctions import basisfunctions
from knot_multiplicity import knot_multiplicity


def Bi(u, a, x, n, n_splines, U, degree):

    i = findspan(n_splines - 1, degree, u, U)
    derivatives(i, u, degree, n, U, x)
    return x[n, a]


def B_prime(u, a, x, n, n_splines, U, degree, band):

    i = findspan(n_splines - 1, degree, u, U)
    derivatives(i, u, degree, n, U, x)

    if (band == 0):
        integrand = x[n, a]**2
    else:
        integrand = x[n, a] * x[n, a + band]
    return integrand


def bvp_boundary_layer(degree, pts):

    # domain
    mink = -1.0
    maxk = 1.0
    # knot vector
    u = np.linspace(mink, maxk, pts)
    # multiplicity of knots according to degree
    U = knot_multiplicity(degree, u)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = 1
    pts = n_splines
    # der[:, k, :] is the kth derivative with the
    # shape of pts x splines
    ders = np.zeros((pts, n + 1, n_splines))

    error = 0.0
    x = np.zeros((n+1, degree + 1))
    size = n_splines - 2

    # boundary layer thickness of O(delta) at x = 1
    delta = 1.0 / 10.0

    # coefficient of y''
    const2 = 1.0
    # coefficient of y'
    const1 = 0.0
    # coefficient of y
    const0 = (1.0 / delta)**2
    # bc at end points
    bc0 = 1.0
    bc1 = 0.0
    # parameters
    n0 = 0
    n1 = 1

    A = np.zeros((size, size))
    b = np.zeros(size)
    f = np.ones(size) * (1.0 / delta)**2

    u = U
    u = np.delete(u, 0)
    u = np.delete(u, -1)

    # position of B_i in respective subintervals
    a = np.arange(degree + 1)[::-1]

    xi = np.zeros(degree + 2)

    # quad rule with B-splines
    for i in range(size):

        # coordinates of subintervals
        for j in range(len(xi)):
            xi[j] = u[i + j]

        for band in range(-degree, degree + 1):
            for j in range(len(a) - abs(band)):

                if (band == 0):
                    # (f_i , B_i)
                    b[i] += f[i] * quad(Bi, xi[j], xi[j+1], args=(a[j], x, n0, n_splines, U, degree))[0]

                if (band >= 0 and (i + band) < size):
                    k = j + band

                    # (B_j', B_i')
                    A[i, i + band] += const2 * quad(B_prime, xi[k], xi[k+1], args=(a[k], x, n1, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    A[i, i + band] += const0 * quad(B_prime, xi[k], xi[k+1], args=(a[k], x, n0, n_splines, U, degree, band))[0]

                elif (band < 0 and (i + band) >= 0):
                    # (B_j', B_i')
                    A[i, i + band] += const2 * quad(B_prime, xi[j], xi[j+1], args=(a[j], x, n1, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    A[i, i + band] += const0 * quad(B_prime, xi[j], xi[j+1], args=(a[j], x, n0, n_splines, U, degree, band))[0]


    # bc fix for first 3 rows
    xi = np.array([u[degree - 1], u[degree]])
    for i in range(degree):
        a = i + 1
        band = -a
        b[i] = b[i] - bc0 * const2 * quad(B_prime, xi[0], xi[1], args=(a, x, n1, n_splines, U, degree, band))[0]
        b[i] = b[i] - bc0 * const0 * quad(B_prime, xi[0], xi[1], args=(a, x, n0, n_splines, U, degree, band))[0]

    # bc fix for last 3 rows
    xi = np.array([u[-(degree + 1)], u[-(degree)]])
    for i in range(1, degree + 1):
        band = i
        a = degree - band
        b[-i] = b[-i] - bc1 * const2 * quad(B_prime, xi[0], xi[1], args=(a, x, n1, n_splines, U, degree, band))[0]
        b[-i] = b[-i] - bc1 * const0 * quad(B_prime, xi[0], xi[1], args=(a, x, n0, n_splines, U, degree, band))[0]

    # solve A x = b
    alpha = la.solve(A, b)

    # fix boundary conditions
    alpha = np.insert(alpha, 0, bc0)
    alpha = np.append(alpha, bc1)

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

    usol = np.linspace(mink, maxk, xx)
    sol = np.ones(xx) - np.sinh((1.0 / delta) * (u + 1.0)) / np.sinh(2.0 * (1.0 / delta))

    error = la.norm(sol - y, np.inf)
    print(error)

#    plt.figure(1, figsize = (12, 10))
#    plt.xticks(fontsize=25)
#    plt.yticks(fontsize=25)
#    plt.plot(u, y, 'o', label='Galerkin')
#    plt.plot(usol, sol, linewidth=2, label='solution')
#    plt.xlim(mink - 0.05, maxk + 0.05)
#    plt.title('Boundary layer problem', fontsize=35)
#    plt.legend(loc='best', fontsize=30)
#    plt.show()

    return error

def bvp_smooth_2(degree, pts):

    # domain
    mink = 0.0
    maxk = 1.0
    # knot vector
    u = np.linspace(mink, maxk, pts)
    # multiplicity of knots according to degree
    U = knot_multiplicity(degree, u)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = 1
    pts = n_splines
    # der[:, k, :] is the kth derivative with the
    # shape of pts x splines
    ders = np.zeros((pts, n + 1, n_splines))

    error = 0.0
    x = np.zeros((n+1, degree + 1))
    size = n_splines - 2

    # coefficient of y''
    const2 = 1.0
    # coefficient of y'
    const1 = 0.0
    # coefficient of y
    const0 = 4.0
    # bc at end points
    bc0 = 0.0
    bc1 = 0.0
    # parameters
    n0 = 0
    n1 = 1

    A = np.zeros((size, size))
    b = np.zeros(size)
    f = -4.0 * np.cosh(np.ones(size))

    u = U
    u = np.delete(u, 0)
    u = np.delete(u, -1)

    # position of B_i in respective subintervals
    a = np.arange(degree + 1)[::-1]

    xi = np.zeros(degree + 2)

    # quad rule with B-splines
    for i in range(size):

        # coordinates of subintervals
        for j in range(len(xi)):
            xi[j] = u[i + j]

        for band in range(-degree, degree + 1):
            for j in range(len(a) - abs(band)):

                if (band == 0):
                    # (f_i , B_i)
                    b[i] += f[i] * quad(Bi, xi[j], xi[j+1], args=(a[j], x, n0, n_splines, U, degree))[0]

                if (band >= 0 and (i + band) < size):
                    k = j + band

                    # (B_j', B_i')
                    A[i, i + band] += const2 * quad(B_prime, xi[k], xi[k+1], args=(a[k], x, n1, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    A[i, i + band] += const0 * quad(B_prime, xi[k], xi[k+1], args=(a[k], x, n0, n_splines, U, degree, band))[0]

                elif (band < 0 and (i + band) >= 0):
                    # (B_j', B_i')
                    A[i, i + band] += const2 * quad(B_prime, xi[j], xi[j+1], args=(a[j], x, n1, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    A[i, i + band] += const0 * quad(B_prime, xi[j], xi[j+1], args=(a[j], x, n0, n_splines, U, degree, band))[0]


    # bc fix for first 3 rows
    xi = np.array([u[degree - 1], u[degree]])
    for i in range(degree):
        a = i + 1
        band = -a
        b[i] = b[i] - bc0 * const2 * quad(B_prime, xi[0], xi[1], args=(a, x, n1, n_splines, U, degree, band))[0]
        b[i] = b[i] - bc0 * const0 * quad(B_prime, xi[0], xi[1], args=(a, x, n0, n_splines, U, degree, band))[0]

    # bc fix for last 3 rows
    xi = np.array([u[-(degree + 1)], u[-(degree)]])
    for i in range(1, degree + 1):
        band = i
        a = degree - band
        b[-i] = b[-i] - bc1 * const2 * quad(B_prime, xi[0], xi[1], args=(a, x, n1, n_splines, U, degree, band))[0]
        b[-i] = b[-i] - bc1 * const0 * quad(B_prime, xi[0], xi[1], args=(a, x, n0, n_splines, U, degree, band))[0]

    # solve A x = b
    alpha = la.solve(A, b)

    # fix boundary conditions
    alpha = np.insert(alpha, 0, bc0)
    alpha = np.append(alpha, bc1)

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

    usol = np.linspace(mink, maxk, xx)
    sol = np.cosh(2.0 * u - 1.0) - np.cosh(1.0)

    error = la.norm(sol - y, np.inf)
    print(error)

#    plt.figure(1, figsize = (12, 10))
#    plt.xticks(fontsize=25)
#    plt.yticks(fontsize=25)
#    plt.plot(u, y, 'o', label='Galerkin')
#    plt.plot(usol, sol, linewidth=2, label='solution')
#    plt.xlim(mink - 0.05, maxk + 0.05)
#    plt.title('Smooth problem 2', fontsize=35)
#    plt.legend(loc='best', fontsize=30)
#    plt.show()

    return error

def bvp_smooth_1(degree, pts):

    # domain
    mink = 0.0
    maxk = 1.0
    # knot vector
    u = np.linspace(mink, maxk, pts)
    # multiplicity of knots according to degree
    U = knot_multiplicity(degree, u)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = 1
    pts = n_splines
    # der[:, k, :] is the kth derivative with the
    # shape of pts x splines
    ders = np.zeros((pts, n + 1, n_splines))

    error = 0.0
    x = np.zeros((n+1, degree + 1))
    size = n_splines - 2

    # coefficient of y''
    const2 = 1.0
    # coefficient of y'
    const1 = 0.0
    # coefficient of y
    const0 = 100.0
    # bc at end points
    bc0 = 1.0
    bc1 = 1.0
    # parameters
    n0 = 0
    n1 = 1

    A = np.zeros((size, size))
    b = np.zeros(size)
    f = np.zeros(size)

    u = U
    u = np.delete(u, 0)
    u = np.delete(u, -1)

    # position of B_i in respective subintervals
    a = np.arange(degree + 1)[::-1]

    xi = np.zeros(degree + 2)

    # quad rule with B-splines
    for i in range(size):

        # coordinates of subintervals
        for j in range(len(xi)):
            xi[j] = u[i + j]

        for band in range(-degree, degree + 1):
            for j in range(len(a) - abs(band)):

                if (band == 0):
                    # (f_i , B_i)
                    b[i] += f[i] * quad(Bi, xi[j], xi[j+1], args=(a[j], x, n0, n_splines, U, degree))[0]

                if (band >= 0 and (i + band) < size):
                    k = j + band

                    # (B_j', B_i')
                    A[i, i + band] += const2 * quad(B_prime, xi[k], xi[k+1], args=(a[k], x, n1, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    A[i, i + band] += const0 * quad(B_prime, xi[k], xi[k+1], args=(a[k], x, n0, n_splines, U, degree, band))[0]

                elif (band < 0 and (i + band) >= 0):
                    # (B_j', B_i')
                    A[i, i + band] += const2 * quad(B_prime, xi[j], xi[j+1], args=(a[j], x, n1, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    A[i, i + band] += const0 * quad(B_prime, xi[j], xi[j+1], args=(a[j], x, n0, n_splines, U, degree, band))[0]


    # bc fix for first 3 rows
    xi = np.array([u[degree - 1], u[degree]])
    for i in range(degree):
        a = i + 1
        band = -a
        b[i] = b[i] - bc0 * const2 * quad(B_prime, xi[0], xi[1], args=(a, x, n1, n_splines, U, degree, band))[0]
        b[i] = b[i] - bc0 * const0 * quad(B_prime, xi[0], xi[1], args=(a, x, n0, n_splines, U, degree, band))[0]

    # bc fix for last 3 rows
    xi = np.array([u[-(degree + 1)], u[-(degree)]])
    for i in range(1, degree + 1):
        band = i
        a = degree - band
        b[-i] = b[-i] - bc1 * const2 * quad(B_prime, xi[0], xi[1], args=(a, x, n1, n_splines, U, degree, band))[0]
        b[-i] = b[-i] - bc1 * const0 * quad(B_prime, xi[0], xi[1], args=(a, x, n0, n_splines, U, degree, band))[0]

    # solve A x = b
    alpha = la.solve(A, b)

    # fix boundary conditions
    alpha = np.insert(alpha, 0, bc0)
    alpha = np.append(alpha, bc1)

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

    usol = np.linspace(mink, maxk, xx)
    sol = np.cosh(10.0 * (usol - 0.5)) / np.cosh(5.0)

    error = la.norm(sol - y, np.inf)
    print(error)

#    plt.figure(1, figsize = (12, 10))
#    plt.xticks(fontsize=25)
#    plt.yticks(fontsize=25)
#    plt.plot(u, y, 'o', label='Galerkin')
#    plt.plot(usol, sol, linewidth=2, label='solution')
#    plt.xlim(mink - 0.05, maxk + 0.05)
#    plt.title('Smooth problem 1', fontsize=35)
#    plt.legend(loc='best', fontsize=30)
#    plt.show()

    return error


def test_bvp_degree_three(degree, problem):
    import csv

    error = []
    pts = np.array([10, 20, 40, 80, 100, 150, 200, 250, 300, 350])
    # pts = np.array([10, 20, 30, 40, 45, 50, 55, 60, 65, 70])

    for i in range(len(pts)):
        if (problem == 1):
            error.append(bvp_smooth_1(degree, pts[i]))
        elif (problem == 2):
            error.append(bvp_smooth_2(degree, pts[i]))
        elif (problem == 3):
            error.append(bvp_boundary_layer(degree, pts[i]))

    temp = zip(pts, error)

    with open('output.csv', "w") as output:
        writer = csv.writer(output, lineterminator='\n', delimiter=',')
        for val in temp:
            writer.writerow(val)


def main():
    test_bvp_degree_three(5, 3)
    #error = bvp_smooth_1(4, 10)
    #error = bvp_smooth_2(4, 10)
    #error = bvp_boundary_layer(4, 10)

if __name__ == "__main__":
    main()
