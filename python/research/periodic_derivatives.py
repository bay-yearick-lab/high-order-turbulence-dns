import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from findspan import findspan
from derivatives import derivatives
from collocation import collocation_points
from wavenumber_periodic import periodic_fix
from knot_multiplicity import periodic_knot_multiplicity


def fourier_periodic_derivatives(pts, oversample):

    # Spatial grid
    m = int(pts * oversample)
    # Width of spatial domain
    L = 2.0 * np.pi
    # Grid spacing
    dx = L / m
    # Grid in physical space
    x = np.arange(-0.5 * m, 0.5 * m) * dx
    # number of derivatives
    deriv = 3

    # Wavenumber grid
    # xi_j = (2 * pi * j) / L,  -m/2 < j < m/2
    xi = np.fft.fftfreq(m) * m * 2 * np.pi / L

    xi_max = np.max(np.abs(xi))
    filtr = np.ones_like(xi)
    if (oversample != 1.0):
        filtr[np.where(np.abs(xi) > int(pts/2 - 1))] = 0.0

    print(filtr)
    # Initial data
    u = np.sin(x)
    uhat = np.fft.fft(u)

    num = np.zeros((deriv+1, len(x)))
    sol = np.zeros((deriv+1, len(x)))

    num[1, :] = np.real(np.fft.ifft(1.0j * xi * uhat * filtr))
    num[2, :] = np.real(np.fft.ifft(- xi**2 * uhat * filtr))
    num[3, :] = np.real(np.fft.ifft(-1.0j * xi**3 * uhat * filtr))

    sol[1, :] = np.cos(x)
    sol[2, :] = -np.sin(x)
    sol[3, :] = -np.cos(x)

    error = sol[1, :] - num[1, :]
    print(abs(error))
    print(la.norm(error, np.inf), np.argmax(abs(error)))

    plt.plot(x, num[1, :],'o')
    plt.plot(x, sol[1, :], linewidth = 2)
    plt.show()


def bspline_periodic_derivatives(degree, pts):
    from collocation_kdv import arbitrary_fix

    # knot vector
    mink = -np.pi
    maxk = np.pi
    U = np.linspace(mink, maxk, pts)
    # knot multiplicity
    U = periodic_knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # collocation pts
    pts = n_splines
    # collocation points
    co = np.zeros((2, pts))
    # number of derivatives
    deriv = 3

    # collocation points
    if (degree%2 == 0):
        pts = n_splines - degree + 2
    else:
        pts = n_splines - degree + 1

    co = np.zeros((2, pts))
    collocation_points(co, U[degree: -(degree)], degree, 'periodic')
    # random u value whom we are trying to find the ith interval it lies in
    pts = pts
    u = co[0, :]
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((pts, deriv + 1, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, deriv, U, ders[j, :, (i - degree):(degree + i)])

    fullD = np.zeros((deriv+1, pts, n_splines))
    D = np.zeros((deriv+1, len(u)-1, len(u)-1))
    num = np.zeros((deriv+1, len(u)))
    sol = np.zeros((deriv+1, len(u)))

    for i in range(deriv+1):
        fullD[i, :, :] = ders[:, i, :]
        D[i, :, :] = periodic_fix(fullD[i, :, :], degree)

    rhs = np.sin(u[1:])
    alpha = la.solve(D[0, :, :], rhs)

    for i in range(1, deriv+1):
        num[i, 1:] = D[i, :, :] @ alpha
        num[i, 0] = num[i, -1]

    sol[1, :] = np.cos(u)
    sol[2, :] = -np.sin(u)
    sol[3, :] = -np.cos(u)

    error = sol[1, :] - num[1, :]
    print(abs(error))
    print(la.norm(error, np.inf), np.argmax(abs(error)))

    plt.plot(u, num[1, :],'o')
    plt.plot(u, sol[1, :], linewidth = 2)
    plt.show()


def periodic_init(mink, maxk, pts, degree, deriv, U2=None):

    U = np.linspace(mink, maxk, pts)
    # knot multiplicity
    U = periodic_knot_multiplicity(degree, U)

    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # collocation pts
    pts = n_splines
    # collocation points
    co = np.zeros((2, pts))

    # collocation points
    if (degree%2 == 0):
        pts = n_splines - degree + 2
    else:
        pts = n_splines - degree + 1

    co = np.zeros((2, pts))
    collocation_points(co, U[degree: -(degree)], degree, 'periodic')
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

    return co[0, :], ders, U


def bspline_dealiasing_periodic_derivatives(degree, pts, oversample):
    from collocation_kdv import arbitrary_fix

    # original points
    og_points = pts
    # oversample
    pts = int(np.round(oversample * og_points))
    # knot vector
    mink = -np.pi
    maxk = np.pi
    # derivatives
    deriv = 1

    u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)

    fullB = ders[:, 0, :].copy()
    fullD = ders[:, 1, :].copy()

    B = periodic_fix(fullB, degree)
    D = periodic_fix(fullD, degree)

    rhs = np.exp(-10.0 * u[1:]**2)
    alpha = la.solve(B, rhs)
    alpha = np.insert(alpha, 0, alpha[-1])

    # interpolation
    xx = og_points

    _, ders2, U2 = periodic_init(mink, maxk, xx, degree, deriv, U)
    interpolation_matrix = periodic_fix(ders2[:, 0, :], degree)
    y = interpolation_matrix @ alpha

    # reducing alpha's by half
    _, ders3, _ = periodic_init(mink, maxk, xx, degree, deriv, U2)
    Bhalf = periodic_fix(ders3[:, 0, :], degree)

    alphahalf = la.solve(Bhalf, y[1:])
    alphahalf = np.insert(alphahalf, 0, alphahalf[-1])

    # extrapolation
    u4, ders4, _ = periodic_init(mink, maxk, pts, degree, deriv, U2)

    extrapolation_matrix = periodic_fix(ders4[:, 0, :], degree, oversample)

    y = extrapolation_matrix @ alphahalf

    sol = np.exp(-10.0 * u4**2)
    error = la.norm(sol - y, np.inf)
    print(error)
    plt.plot(u4, y,'o')
    plt.plot(u4, sol, linewidth = 2)
    plt.show()


def main():
    bspline_dealiasing_periodic_derivatives(5, 400, 2)
    # bspline_periodic_derivatives(5, 20)
    # fourier_periodic_derivatives(20, 1)


if __name__ == "__main__":
    main()
