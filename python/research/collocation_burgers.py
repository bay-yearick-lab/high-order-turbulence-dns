import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import matplotlib.pyplot as plt
from findspan import findspan
from derivatives import derivatives
from collocation import collocation_points
from knot_multiplicity import knot_multiplicity


def runge_kutta(stage):
    if (stage == 6):
        return np.array([1.0/9.0, 2.0/5.0, 0.0, 4.0/9.0, 0.0, 2.0/45.0])
    elif (stage == 8):
        return np.array([2.0/15.0, 2.0/7.0, 2.0/9.0, 4.0/15.0, 0.0,
                        4.0/45.0, 0.0, 1.0/315.0])
    elif (stage == 10):
        return np.array([71./525, 22./81, 4./15, 4./21, 2./27, 4./75,
                        0.0, 8./945, 0.0, 2./14175])


def htan(x, nu, t):
    f1 = 1.0
    f2 = 0.0
    f = f1 - f2
    c = 0.5 * (f1 + f2)
    return c - 0.5 * f * np.tanh(0.25/nu * f * (x - c * t))


def f(x, nu, t):
    return htan(x, nu, t)


def construct_M(diffusion, convection, u):

    convection = - np.diag(u) @ convection
    M = convection + diffusion


    firstcolumn = M[:, 0]
    firstcolumn = np.delete(firstcolumn, 0)
    firstcolumn = np.delete(firstcolumn, -1)
    lastcolumn = M[:, -1]
    lastcolumn = np.delete(lastcolumn, 0)
    lastcolumn = np.delete(lastcolumn, -1)

    M = np.delete(M, 0, 0)
    M = np.delete(M, -1, 0)
    M = np.delete(M, 0, 1)
    M = np.delete(M, -1, 1)

    return M, firstcolumn, lastcolumn


def collocation_burgers(degree, pts):

    # knot vector
    mink = -4.0
    maxk = 4.0
    U = np.linspace(mink, maxk, pts)
    # knot multiplicity
    U = knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # collocation pts
    pts = n_splines
    # collocation points
    co = np.zeros((2, pts))
    # number of derivatives
    n = 2
    # time
    t = 0.0
    # viscosity
    nu = 0.1
    # runge kutta stage
    stage = 6

    collocation_points(co, U, degree)

    u = co[0, :]
    ders = np.zeros((pts, n + 1, n_splines))
    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    fullB = ders[:, 0, :]
    fullBinv = la.inv(fullB)
    fullD = ders[:, 1, :]
    fullL = ders[:, 2, :]

    fullu_n = np.zeros((stage, len(u)))
    u_n = np.zeros((stage, len(u) - 2))
    fullu_n[0, :] = f(u, nu, t)
    u_n[0, :] = fullu_n[0, 1:-1].copy()

    rk = runge_kutta(stage)

    diffusion = nu * (fullL @ fullBinv)
    convection = fullD @ fullBinv
    M, firstcolumn, lastcolumn = construct_M(diffusion, convection, fullu_n[0, :])

#    eigvals = la.eigvals(M)
#    print('max real = %f'%max(np.real(eigvals)))
#    print('min real = %f'%min(np.real(eigvals)))
#    print('dt = %f'% (-2.0 / min(np.real(eigvals))))
#    plt.plot(np.real(eigvals), np.imag(eigvals), 'o')
#    plt.show()
#    quit()

    # time
    t = 0.0
    # final time
    T = 1.0
    # time-step
    dt = 5.0e-4

    print('dt = %f' % dt)

    sol = f(u, nu, t)

#    xx = np.linspace(mink, maxk, 5000)
#    plt.ion()
#    plt.figure(1, figsize=(18, 15))
#    plt.hold(True)
#    udata, = plt.plot(u, f(u, nu, t), 'o', markersize=10, label='Numerical')
#    uline, = plt.plot(u, f(u, nu, t), '-', linewidth=2)
#    exact, = plt.plot(xx, f(xx, nu, t), '-g', linewidth=4, zorder=-10, label='Exact')
#    error, = plt.plot(
#             u, np.zeros_like(u), '-r', linewidth=2, zorder=-10,
#             label='$L_\infty$ error = %0.2E' % la.norm(sol-fullu_n[0, :], np.inf))
#    plt.legend()
#    plt.pause(0.001)

    N = int(T / dt)

    for time in range(N):

        # runge kutta SSP
        for m in range(stage):
            t += (rk[m] * dt)

            bc0 = f(u[0], nu, t)
            bc1 = f(u[-1], nu, t)

            M, firstcolumn, lastcolumn = construct_M(diffusion, convection, fullu_n[m, :])

            firstcolumn = 0.5 * dt * firstcolumn * bc0
            lastcolumn = 0.5 * dt * lastcolumn * bc1

            u_n[m, :] = np.delete(np.delete(fullu_n[m, :], 0), -1)

            if ((m+1) != stage):
                u_n[m + 1, :] = u_n[m, :] + 0.5 * dt * (M @ u_n[m, :]) + firstcolumn + lastcolumn

                a = np.insert(u_n[m + 1, :], 0, bc0)
                fullu_n[m + 1, :] = np.append(a, bc1)
            else:
                u_n[0, :] = np.dot(rk[:m-1], u_n[:m-1, :]) + \
                            rk[-1] * (u_n[m, :] + 0.5 * dt * (M @ u_n[m, :]) + firstcolumn + lastcolumn)

                a = np.insert(u_n[0, :], 0, bc0)
                fullu_n[0, :] = np.append(a, bc1)

        sol = f(u, nu, t)
        Linf = la.norm(sol - fullu_n[0, :], np.inf)
        assert(Linf < 1e2)

#        print('Time = %f' % t)
#        udata.set_ydata(fullu_n[0, :])
#        uline.set_ydata(fullu_n[0, :])
#        exact.set_ydata(f(xx, nu, t))
#        error.set_ydata(sol - fullu_n[0, :])
#        error.set_label('$L_\infty$ error = %0.2E' % Linf)
#
#        plt.ylabel('$u(x)$', fontsize=55)
#        plt.xlabel('$x$', fontsize=55)
#        plt.xticks(fontsize=30)
#        plt.yticks(fontsize=30)
#        plt.title('%dth-Order B-Spline \n (Burgers\' Equation, T = %0.2f)'
#                  % (degree, t), fontsize=45)
#        plt.legend(fontsize=40, loc='best')
#        plt.ylim(-0.2 + min(fullu_n[0, :]), 0.2 + max(fullu_n[0, :]))
#        plt.pause(0.001)
#        plt.draw()
#        if (t > T):
#            plt.ioff()
#            plt.show()

    a = np.argmax(abs(sol - fullu_n[0, :]))
    print(pts, Linf, a)
    return Linf


def test_collocation_burgers(degree):
    import csv

    error = []
    pts = np.array([20, 40, 60, 80, 100, 120, 140, 160, 180])

    for i in range(len(pts)):
        error.append(collocation_burgers(degree, pts[i]))

    temp = zip(pts, error)

    with open('collocation-burgers-degree%d.csv'% degree, "w") as output:
        writer = csv.writer(output, lineterminator='\n', delimiter=',')
        for val in temp:
            writer.writerow(val)


def main():
    # collocation_burgers(5, 180)
    test_collocation_burgers(5)


if __name__ == "__main__":
    main()
