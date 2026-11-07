import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
from findspan import findspan
from derivatives import derivatives
from collocation import collocation_points
from wavenumber_periodic import periodic_fix
from knot_multiplicity import periodic_knot_multiplicity
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def objective(lamda, B, I, P, u):
    A = B @ la.inv(B.T @ B + lamda * P) @ B.T
    V = la.norm(u - A @ u, 2)**2
    V = V / np.trace(I - A)**2
    return V

def objective_vector(lamda, B, I, P, u):
    A = B @ la.inv(B.T @ B + np.diag(lamda) @ P) @ B.T
    V = la.norm(u - A @ u, 2)**2
    V = V / np.trace(I - A)**2
    return V


def runge_kutta(stage):
    if (stage == 6):
        return np.array([1.0/9.0, 2.0/5.0, 0.0, 4.0/9.0, 0.0, 2.0/45.0])
    elif (stage == 8):
        return np.array([2.0/15.0, 2.0/7.0, 2.0/9.0, 4.0/15.0, 0.0,
                        4.0/45.0, 0.0, 1.0/315.0])
    elif (stage == 10):
        return np.array([71./525, 22./81, 4./15, 4./21, 2./27, 4./75,
                        0.0, 8./945, 0.0, 2./14175])


def f(x, t):
    A = 25.0; B = 16.0
    u = 3.0 * A**2 / np.cosh(0.5 * (A * (x + 2.0)))**2 \
      + 3.0 * B**2 / np.cosh(0.5 * (B * (x + 1.0)))**2
    return u


def f1(x):
    return np.exp(-100.0 * (x)**2)


def sech_pulse(x):
    return 1.0 / np.cosh(4 * (x - np.mean(x)))


def compact_filter(pts, alpha):
    Erow = np.zeros(pts)
    Irow = np.zeros(pts)
    # 8th order compact
    a = 1/128.0 * (93.0 + 70 * alpha)
    b = 1.0/16.0 * (7.0 + 18.0 * alpha)
    c = 1.0/32.0 * (-7.0 + 14.0 *alpha)
    d = 1.0/16.0 * (1.0 - 2.0 * alpha)
    e = 1.0/128.0 * (-1.0 + 2.0 *alpha)

    Erow[-4] = e/2; Erow[-3] = d/2; Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a
    Erow[4] = e/2; Erow[3] = d/2; Erow[2] = c/2; Erow[1] = b/2
    explicit8 = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit8 = np.zeros((pts, pts))

    for i in range(pts):
        explicit8[i, :] = np.roll(Erow, i)
        implicit8[i, :] = np.roll(Irow, i)
    return explicit8, implicit8


def ic_burgers_turbulence(pts):
    psi = np.random.uniform(0.0, 1.0, pts)
    k = np.linspace(-np.pi*pts/(2.0*np.pi), np.pi*pts/(2.0*np.pi), pts)
    A = 2 * 10.0**(-5) / (3.0 * np.sqrt(np.pi))
    E = A * k**4 * np.exp(-(k/10.0)**2)

    median = int(pts/2)
    psi[:median] = -psi[median:]

    uhat = np.sqrt(2 * E) * np.exp(1j * 2.0 * np.pi * psi)

    u = np.real(np.fft.fftshift(np.fft.ifft(uhat)))
    return u


def construct_M(dispersion, convection, u):
     return - np.diag(u) @ convection - dispersion


def plot_frame(i, uline, axes, x, frames, tt, pts):
    uline.set_data(x, frames[i])
    axes.set_title('%d points, Time = %e'%(pts, tt[i]), fontsize=40)
    axes.set_xlim((min(x), max(x)))
    #axes.set_ylim((-100,3000))

    axes.set_ylabel('$u(x, t)$', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)


def arbitrary_fix(M, degree):

    if (degree%2 == 0):
        M = np.delete(M, -1, 0)
        M = np.delete(M, 0, 0)
        m, n = M.shape
        row = np.zeros((degree - 1, m))

        for i in range(degree):
            M = np.delete(M, -1, 1)
            if (i < (degree-1)):
                n = degree - (i + 1)
                row[i, :] = M[-(i+2), n:]
                row[i, :] = np.roll(row[i, :], n)

        for i in range(1, degree):
            M[-i, :] = row[i-1, :]

    else:
        M = np.delete(M, -1, 0)
        m, n = M.shape
        row = np.zeros((degree, m))
        for i in range(degree):
            M = np.delete(M, -1, 1)
            if (i < (degree)):
                n = degree - (i + 1)
                row[i, :] = M[-(i + 2), n:]
                row[i, :] = np.roll(row[i, :], n+1)

        for i in range(1, degree+1):
            M[-i, :] = row[i-1, :]
    return M


def collocation_burgers_smooth_vector(degree, pts, p_deriv, convective_type):
    from periodic_derivatives import periodic_init
    from exact_periodic_burgers import exact_burgers
    from smooth_spline import smooth_spline_matrix_init
    from smooth_spline import smooth_spline_rhs_init
    from scipy.optimize import minimize_scalar
    from scipy.optimize import minimize

    # save original number of points
    og_points = pts
    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2
    # time
    t = 0.0
    # runge kutta stage
    stage = 6
    # time
    t = 0.0
    # final time
    T = 4.0
    # time-step
    dt = 1.0e-3
    # diffusion
    nu = 1.0 / 80.0

    u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)

    B = periodic_fix(ders[:, 0, :], degree)
    convection = periodic_fix(ders[:, 1, :], degree)
    diffusion = nu * periodic_fix(ders[:, 2, :], degree)

    u_n = np.zeros((stage, len(u)))
    # Initial data
    u_n[0, :] = sech_pulse(u)

    rk = runge_kutta(stage)

    print('%s, %d degree, %d points'%(convective_type, degree, pts))

    lamda = 0.4e-6
    A, firstcolumn, lastcolumn, fullB, P = \
            smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, p_deriv, 'periodic')

    BT = fullB.T

    N = int(round(T / dt))

    LU = sla.lu_factor(B)
    smooth_LU = sla.lu_factor(A)

    plt.rc('text', usetex=True)
    fig = plt.figure(1, figsize=(9, 4))
    ax = fig.add_subplot(111)
    uline, = ax.plot(u, u_n[0, :], '-', linewidth=2)

    penalty = True

    shape, _ = BT.shape
    I = np.eye(shape)

    Prow = P[degree, :]

    for i in range(shape):
        P[i, :] = np.roll(Prow, -degree + i)

    lamda = np.zeros(shape)
    for time in range(N):

        lamda = 2.0e-7 + 0.0 * lamda
        #res = minimize(objective_vector, lamda, args=(fullB, I, P, u_n[0, :]), \
        #        method='L-BFGS-B')

        #if res.success != True:
        #    print("Optimization failed")
        #    lamda = 5.0e-7 * lamda
        #else:
        #    lamda = res.x

        #lamda[lamda<0.0] = 5.0e-7
        #lamda[lamda>1.0] = 5.0e-7

        umax = 0.5 * sum(abs(u_n[0, :]))**2
        # runge kutta SSP
        for m in range(stage):
            t += (rk[m] * dt)

            # smoothing spline
            smooth_LU = sla.lu_factor((BT @ fullB + np.diag(lamda) @ P)[1:, 1:])
            rhs, bc = smooth_spline_rhs_init(u_n[m, :], firstcolumn, lastcolumn, BT, 'periodic')
            alpha = sla.lu_solve(smooth_LU, rhs)
            u_n[m, 1:] = B @ alpha

            if (convective_type == 'Advective'):
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                                    + 0.5 * dt * (diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                              + 0.5 * dt * (diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Divergence'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * (
                                    - 0.5 * (convection @ alpha2) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              - 0.5 * (convection @ alpha2) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Skew-symmetric'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * ( \
                                    -0.25 * (convection @ alpha2) \
                                    -0.50 * (u_n[m, 1:] * (convection @ alpha)) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              -0.25 * (convection @ alpha2) \
                              -0.50 * (u_n[0, 1:] * (convection @ alpha)) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]


        #umax2 = 0.5 * sum(abs(u_n[0, :]))**2
        #if umax2 > umax:
        #    penalty = True
        #    print('Time = %e, L_inf = %e' % (t, umax2))
        #else:
        #    penalty = False

        print('Time = %e, Norm = %e' % (t, la.norm(u_n[0, :], np.inf)))
        assert(la.norm(u_n[0, :], np.inf) < 1e2)

#        # saving frames for video
#        if np.mod(time, nplt) == 0:
#            frames.append(u_n[0, :].copy())
#            tt.append(t)

        uline.set_data(u, u_n[0, :])
        ax.set_ylabel('$u(x, t)$', fontsize=30)
        ax.set_xlabel('$x$', fontsize=30)
        plt.ylim(top=0.65)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        ax.set_title('B-splines (with filtering), points = 32', fontsize=40)

        if (time == N-1):
            plt.tight_layout()
#            plt.savefig('op-burgers.eps', transparent=True, dpi=300)
            plt.show()
#            t = np.linspace(0.0, T, len(lam_all))
#            plt.figure(2, figsize=(9, 4))
#            plt.ylabel('$\\lambda$', fontsize=30)
#            plt.xlabel('$t$', fontsize=30)
#            plt.xticks(fontsize=20)
#            plt.yticks(fontsize=20)
#            plt.plot(t, lam_all, '-', linewidth=2)
#            plt.title('Value of $\\lambda$ vs $t$', fontsize=40)
#            plt.tight_layout()
#            plt.savefig('op-lamda-burgers.eps', transparent=True, dpi=300)
#            quit()


#    # movie
#    im_ani = matplotlib.animation.FuncAnimation(fig, plot_frame, frames=len(frames),
#            fargs=(uline, axes, u, frames, tt, pts), interval=20)
#    im_ani.save('burgers.mp4', writer=writer)

#    exact = exact_burgers(u, og_points, R, T)
#    Linf = la.norm(exact - u_n[0, :], np.inf)
#    L2 = np.sqrt(sum(abs(exact - u_n[0, :])**2) / len(exact))
#
#    print(convective_type, '&  %0.2f'%Linf, '&  %0.2f'%L2)

    data = np.array([u, u_n[0, :]]).T
    with open("scalar.out", 'w+') as file:
            np.savetxt(file, data, fmt=['%20.10f','%19.10f'])
    plt.show()

    return u, u_n[0, :]


def collocation_burgers(degree, pts, p_deriv, convective_type):
    from periodic_derivatives import periodic_init
    from exact_periodic_burgers import exact_burgers
    from smooth_spline import smooth_spline_matrix_init
    from smooth_spline import smooth_spline_rhs_init
    from scipy.optimize import minimize_scalar
    from scipy.optimize import minimize

    # Set up formatting for the movie files
#    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='burgers'), bitrate=1800)

    # save original number of points
    og_points = pts
    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2
    # time
    t = 0.0
    # runge kutta stage
    stage = 6
    # time
    t = 0.0
    # final time
    T = 4.0
    # time-step
    dt = 1.0e-3
    # diffusion
    nu = 1.0 / 80.0

    u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)

    B = periodic_fix(ders[:, 0, :], degree)
    convection = periodic_fix(ders[:, 1, :], degree)
    diffusion = nu * periodic_fix(ders[:, 2, :], degree)

#    Binv = la.inv(B)
#    M = -convection @ Binv
#
#    eigvals = la.eigvals(M)
#    plt.plot(np.real(eigvals), np.imag(eigvals), 'o')
#    plt.xlabel('$Re$', fontsize=35)
#    plt.ylabel('$Im$', fontsize=35)
#    plt.axis('equal')
#    plt.show()
#    print(max(np.real(eigvals)))
#    print(min(np.real(eigvals)))
#    quit()

    u_n = np.zeros((stage, len(u)))
    # Initial data
    u_n[0, :] = sech_pulse(u)

    rk = runge_kutta(stage)

    print('%s, %d degree, %d points'%(convective_type, degree, pts))

    lamda = 1e-6
    A, firstcolumn, lastcolumn, fullB, P = \
            smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, p_deriv, 'periodic')

    BT = fullB.T

#    # movie
#    plt.rc('text', usetex=True)
#    fig = plt.figure(figsize=(9, 4))
#    axes = fig.add_subplot(211)
#    uline, = axes.plot(u, u_n[0, :], '-', linewidth=2)
#    frames = [u_n[0, :].copy()]
#    tt = [0]
#    num_plots = 200
#    nplt = np.floor((T / num_plots) / dt)

    N = int(round(T / dt))

    LU = sla.lu_factor(B)
    smooth_LU = sla.lu_factor(A)

    plt.rc('text', usetex=True)
    fig = plt.figure(1, figsize=(9, 4))
    ax = fig.add_subplot(111)
    uline, = ax.plot(u, u_n[0, :], '-', linewidth=2)

    penalty = False

    shape, _ = BT.shape
    I = np.eye(shape)

    Prow = P[degree, :]

    for i in range(shape):
        P[i, :] = np.roll(Prow, -degree + i)

    lam_all = np.zeros(N)
    for time in range(N):

        if time == N-1:
            penalty = True
            lamda=1e-5
            #lamda = minimize(objective, ((1.0e-5),), args=(fullB, I, P, u_n[0, :]), \
            #        method='BFGS').x
            #if (lamda <0.0 or lamda>1.0):
            #    lamda = 0.0

            #lam_all[time] = lamda

        umax = 0.5 * sum(abs(u_n[0, :]))**2
        # runge kutta SSP
        for m in range(stage):
            t += (rk[m] * dt)

            if (np.isnan(u_n[m, :]).any() or np.isinf(u_n[m, :]).any()):
                KE_time = np.linspace(0.0, time * dt, len(KE))
                return KE_time, KE

            if penalty == False:
                alpha = sla.lu_solve(LU, u_n[m, 1:])
            else:
                # smoothing spline
                smooth_LU = sla.lu_factor((BT @ fullB + lamda * P)[1:, 1:])
                rhs, bc = smooth_spline_rhs_init(u_n[m, :], firstcolumn, lastcolumn, BT, 'periodic')
                alpha = sla.lu_solve(smooth_LU, rhs)
                u_n[m, 1:] = B @ alpha

            if (convective_type == 'Advective'):
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                                    + 0.5 * dt * (diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                              + 0.5 * dt * (diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Divergence'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * (
                                    - 0.5 * (convection @ alpha2) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              - 0.5 * (convection @ alpha2) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Skew-symmetric'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * ( \
                                    -0.25 * (convection @ alpha2) \
                                    -0.50 * (u_n[m, 1:] * (convection @ alpha)) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              -0.25 * (convection @ alpha2) \
                              -0.50 * (u_n[0, 1:] * (convection @ alpha)) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

        #umax2 = 0.5 * sum(abs(u_n[0, :]))**2
        #if umax2 > umax:
        #    penalty = True
        #    print('Time = %e, L_inf = %e' % (t, umax2))
        #else:
        #    penalty = False

        print('Time = %e, Norm = %e' % (t, la.norm(u_n[0, :], np.inf)))
        assert(la.norm(u_n[0, :], np.inf) < 1e2)

#        # saving frames for video
#        if np.mod(time, nplt) == 0:
#            frames.append(u_n[0, :].copy())
#            tt.append(t)

        uline.set_data(u, u_n[0, :])
        ax.set_ylabel('$u(x, t)$', fontsize=30)
        ax.set_xlabel('$x$', fontsize=30)
        plt.ylim(top=0.65)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        ax.set_title('B-splines, points = 256', fontsize=40)

        if (time == N-1):
            plt.tight_layout()
#            plt.savefig('op-burgers.eps', transparent=True, dpi=300)

#            t = np.linspace(0.0, T, len(lam_all))
#            plt.figure(2, figsize=(9, 4))
#            plt.ylabel('$\\lambda$', fontsize=30)
#            plt.xlabel('$t$', fontsize=30)
#            plt.xticks(fontsize=20)
#            plt.yticks(fontsize=20)
#            plt.plot(t, lam_all, '-', linewidth=2)
#            plt.title('Value of $\\lambda$ vs $t$', fontsize=40)
#            plt.tight_layout()
#            plt.savefig('op-lamda-burgers.eps', transparent=True, dpi=300)
#            quit()

    data = np.array([u, u_n[0, :]]).T
    with open("burgers.out", 'w+') as file:
            np.savetxt(file, data, fmt=['%20.10f','%19.10f'])
    plt.show()

#    # movie
#    im_ani = matplotlib.animation.FuncAnimation(fig, plot_frame, frames=len(frames),
#            fargs=(uline, axes, u, frames, tt, pts), interval=20)
#    im_ani.save('burgers.mp4', writer=writer)

#    exact = exact_burgers(u, og_points, R, T)
#    Linf = la.norm(exact - u_n[0, :], np.inf)
#    L2 = np.sqrt(sum(abs(exact - u_n[0, :])**2) / len(exact))
#
#    print(convective_type, '&  %0.2f'%Linf, '&  %0.2f'%L2)

    return u, u_n[0, :]



def collocation_burgers_compact(degree, pts, p_deriv, convective_type):
    from periodic_derivatives import periodic_init
    from exact_periodic_burgers import exact_burgers
    from smooth_spline import smooth_spline_matrix_init
    from smooth_spline import smooth_spline_rhs_init
    from scipy.optimize import minimize_scalar
    from scipy.optimize import minimize

    # Set up formatting for the movie files
#    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='burgers'), bitrate=1800)

    # save original number of points
    og_points = pts
    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2
    # time
    t = 0.0
    # runge kutta stage
    stage = 6
    # time
    t = 0.0
    # final time
    T = 4.0
    # time-step
    dt = 1.0e-3
    # diffusion
    nu = 1.0 / 80.0

    u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)

    B = periodic_fix(ders[:, 0, :], degree)
    convection = periodic_fix(ders[:, 1, :], degree)
    diffusion = nu * periodic_fix(ders[:, 2, :], degree)

#    Binv = la.inv(B)
#    M = -convection @ Binv
#
#    eigvals = la.eigvals(M)
#    plt.plot(np.real(eigvals), np.imag(eigvals), 'o')
#    plt.xlabel('$Re$', fontsize=35)
#    plt.ylabel('$Im$', fontsize=35)
#    plt.axis('equal')
#    plt.show()
#    print(max(np.real(eigvals)))
#    print(min(np.real(eigvals)))
#    quit()

    u_n = np.zeros((stage, len(u)))
    # Initial data
    u_n[0, :] = sech_pulse(u)

    rk = runge_kutta(stage)

    print('%s, %d degree, %d points'%(convective_type, degree, pts))

    lamda = 1e-6
    A, firstcolumn, lastcolumn, fullB, P = \
            smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, p_deriv, 'periodic')

    BT = fullB.T

#    # movie
#    plt.rc('text', usetex=True)
#    fig = plt.figure(figsize=(9, 4))
#    axes = fig.add_subplot(211)
#    uline, = axes.plot(u, u_n[0, :], '-', linewidth=2)
#    frames = [u_n[0, :].copy()]
#    tt = [0]
#    num_plots = 200
#    nplt = np.floor((T / num_plots) / dt)

    N = int(round(T / dt))

    LU = sla.lu_factor(B)
    smooth_LU = sla.lu_factor(A)

    plt.rc('text', usetex=True)
    fig = plt.figure(1, figsize=(9, 4))
    ax = fig.add_subplot(111)
    uline, = ax.plot(u, u_n[0, :], '-', linewidth=2)

    penalty = False

    shape, _ = BT.shape
    I = np.eye(shape)

    Prow = P[degree, :]

    for i in range(shape):
        P[i, :] = np.roll(Prow, -degree + i)

    lam_all = np.zeros(N)

    explicit, implicit = compact_filter(len(u_n[0, :]), 0.4995)
    LU_filter = sla.lu_factor(implicit)
    og_ke = np.zeros(N)
    for time in range(N):

        if time == N-1:
            penalty = True
            lamda=1e-5
            #lamda = minimize(objective, ((1.0e-5),), args=(fullB, I, P, u_n[0, :]), \
            #        method='BFGS').x
            #if (lamda <0.0 or lamda>1.0):
            #    lamda = 0.0

            #lam_all[time] = lamda

        umax = 0.5 * sum(abs(u_n[0, :]))**2
        # runge kutta SSP
        for m in range(stage):
            t += (rk[m] * dt)

            if (np.isnan(u_n[m, :]).any() or np.isinf(u_n[m, :]).any()):
                KE_time = np.linspace(0.0, time * dt, len(KE))
                return KE_time, KE

            if penalty == False:
                alpha = sla.lu_solve(LU, u_n[m, 1:])
            else:
                # smoothing spline
                smooth_LU = sla.lu_factor((BT @ fullB + lamda * P)[1:, 1:])
                rhs, bc = smooth_spline_rhs_init(u_n[m, :], firstcolumn, lastcolumn, BT, 'periodic')
                alpha = sla.lu_solve(smooth_LU, rhs)
                u_n[m, 1:] = B @ alpha

            if (convective_type == 'Advective'):
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                                    + 0.5 * dt * (diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                              + 0.5 * dt * (diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Divergence'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * (
                                    - 0.5 * (convection @ alpha2) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              - 0.5 * (convection @ alpha2) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Skew-symmetric'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * ( \
                                    -0.25 * (convection @ alpha2) \
                                    -0.50 * (u_n[m, 1:] * (convection @ alpha)) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              -0.25 * (convection @ alpha2) \
                              -0.50 * (u_n[0, 1:] * (convection @ alpha)) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

        og_ke[time] = 0.5 * sum(u_n[0, :]**2)
        u_n[0, :] = sla.lu_solve(LU_filter, explicit @ u_n[0, :])
        og_ke[time] = (og_ke[time] - 0.5 * sum(u_n[0, :]**2))

        #umax2 = 0.5 * sum(abs(u_n[0, :]))**2
        #if umax2 > umax:
        #    penalty = True
        #    print('Time = %e, L_inf = %e' % (t, umax2))
        #else:
        #    penalty = False

        print('Time = %e, Norm = %e' % (t, la.norm(u_n[0, :], np.inf)))
        assert(la.norm(u_n[0, :], np.inf) < 1e2)

#        # saving frames for video
#        if np.mod(time, nplt) == 0:
#            frames.append(u_n[0, :].copy())
#            tt.append(t)

        uline.set_data(u, u_n[0, :])
        ax.set_ylabel('$u(x, t)$', fontsize=30)
        ax.set_xlabel('$x$', fontsize=30)
        plt.ylim(top=0.65)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        ax.set_title('B-splines, points = 32', fontsize=40)

        if (time == N-1):
            plt.tight_layout()
#            plt.savefig('op-burgers.eps', transparent=True, dpi=300)

#            t = np.linspace(0.0, T, len(lam_all))
#            plt.figure(2, figsize=(9, 4))
#            plt.ylabel('$\\lambda$', fontsize=30)
#            plt.xlabel('$t$', fontsize=30)
#            plt.xticks(fontsize=20)
#            plt.yticks(fontsize=20)
#            plt.plot(t, lam_all, '-', linewidth=2)
#            plt.title('Value of $\\lambda$ vs $t$', fontsize=40)
#            plt.tight_layout()
#            plt.savefig('op-lamda-burgers.eps', transparent=True, dpi=300)
#            quit()

#    plt.semilogy(np.linspace(0, T, N), og_ke)
#    data = np.array([np.linspace(0, T, N), og_ke]).T
#    with open("compact_burgers.out", 'w+') as file:
#            np.savetxt(file, data, fmt=['%20.10f','%19.10f'])
    data = np.array([u, u_n[0, :]]).T
    with open("constant2.out", 'w+') as file:
            np.savetxt(file, data, fmt=['%20.10f','%19.10f'])
    plt.show()

#    # movie
#    im_ani = matplotlib.animation.FuncAnimation(fig, plot_frame, frames=len(frames),
#            fargs=(uline, axes, u, frames, tt, pts), interval=20)
#    im_ani.save('burgers.mp4', writer=writer)

#    exact = exact_burgers(u, og_points, R, T)
#    Linf = la.norm(exact - u_n[0, :], np.inf)
#    L2 = np.sqrt(sum(abs(exact - u_n[0, :])**2) / len(exact))
#
#    print(convective_type, '&  %0.2f'%Linf, '&  %0.2f'%L2)

    return u, u_n[0, :]

def collocation_kdv(degree, pts):

    # Set up formatting for the movie files
    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='KdV'), bitrate=1800)

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
    n = 3
    # time
    t = 0.0
    # runge kutta stage
    stage = 6
    # time
    t = 0.0
    # final time
    T = 6.0e-3
    # time-step
    dt = 1.0e-7

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
    ders = np.zeros((pts, n + 1, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    fullB = ders[:, 0, :]
    fullD = ders[:, 1, :]
    fullL = ders[:, 3, :]

    L = arbitrary_fix(fullL, degree)
    B = periodic_fix(fullB, degree)
    D = periodic_fix(fullD, degree)

    convection = D
    dispersion = L

    u_n = np.zeros((stage, len(u)-1))
    u_n[0, :] = f(u[1:], t)
    rk = runge_kutta(stage)

#    M = construct_M(dispersion, convection, u_n[0, :])
#
#    eigvals = la.eigvals(M)
#    plt.plot(np.real(eigvals), np.imag(eigvals), 'o')
#    plt.axis('equal')
#    plt.show()
#    print(max(np.real(eigvals)))
#    print('dt = %e'%(2.0 / max(abs(np.imag(eigvals)))))
#    quit()

#    plt.rc('text', usetex=True)
#    fig = plt.figure(1, figsize=(9, 4))
#    ax = fig.add_subplot(111)
#    uline, = ax.plot(u[1:], u_n[0, :], '-', linewidth=2)

    print('%d points'%pts)

    # movie
    plt.rc('text', usetex=True)
    fig = plt.figure(figsize=(9, 4))
    axes = fig.add_subplot(211)
    uline, = axes.plot(u[1:], u_n[0, :], '-', linewidth=2)

    frames = [u_n[0, :].copy()]
    tt = [0]
    num_plots = 50
    nplt = np.floor((T / num_plots) / dt)

    N = int(round(T / dt))

    LU = sla.lu_factor(B)

    diff = np.diff(u)
    omega = np.diag(diff[0] * np.ones(len(u_n[0, :])))
    KE0 = 0.5 * np.inner(u_n[0, :], (omega @ u_n[0, :]))
    KE = []

    for time in range(N):

        # runge kutta SSP
        for m in range(stage):
            t += (rk[m] * dt)

            if (np.isnan(u_n[m, :]).any() or np.isinf(u_n[m, :]).any()):
                KE_time = np.linspace(0.0, time * dt, len(KE))
                return KE_time, KE

            alpha = sla.lu_solve(LU, u_n[m, :])

            if ((m+1) != stage):
                u_n[m + 1, :] = u_n[m, :]  \
                              - 0.5 * dt * u_n[m, :] * (convection @ alpha) \
                              - 0.5 * dt * (dispersion @ alpha)
            else:
                u_n[0, :] = np.dot(rk[:m-1], u_n[:m-1, :]) + \
                            rk[-1] * (u_n[m, :] \
                          - 0.5 * dt * u_n[m, :] * (convection @ alpha) \
                          - 0.5 * dt * (dispersion @ alpha))

#        print('Time = %e' % t)

        KE.append(abs(0.5 * np.inner(u_n[0, :], (omega @ u_n[0, :])) - KE0))

        # saving frames for video
        if np.mod(time, nplt) == 0:
            frames.append(u_n[0, :].copy())
            tt.append(t)

#        print(u_n[0, :])
#        uline.set_data(u[1:], u_n[0, :])
#        ax.set_ylabel('$u(x, t)$', fontsize=30)
#        plt.xticks(fontsize=20)
#        plt.yticks(fontsize=20)
#        ax.set_title('%d points, B-spline degree %d'%(pts, degree), fontsize=40)
#
#        if (time == int(4.0e-3/dt)):
#            plt.savefig('degree%d-b-spline-%d.eps'%(degree, pts), transparent=True, dpi=300)
#            quit()
#
#    fig2 = plt.figure(figsize=(9, 4))
#    plt.ylabel('$|KE(t) - KE_0|$', fontsize=30)
#    plt.xlabel('$t$', fontsize=30)
#    plt.xticks(fontsize=20)
#    plt.yticks(fontsize=20)
#    plt.title('%d points (unfiltered)' % pts, fontsize=40)
#    plt.plot(np.linspace(0.0, T, N), KE)
#    plt.savefig('KE-degree%d-b-spline-%d.eps'%(degree, pts), transparent=True, dpi=300)

    # movie
    im_ani = matplotlib.animation.FuncAnimation(fig, plot_frame, frames=len(frames),
            fargs=(uline, axes, u[1:], frames, tt, pts), interval=20)
    im_ani.save('kdv.mp4', writer=writer)

    KE_time = np.linspace(0.0, T, N)

    return KE_time, KE


def discrete_KE(degree):
    import csv

    pts = np.array([32])
    string = 'unfilter'
    convection_type = ['Skew-symmetric']
    for i in range(len(convection_type)):
        a, b  = collocation_burgers(degree, pts[0], convection_type[i])

#       plt.rc('text', usetex=True)
#       plt.figure(3, figsize=(18, 16))
#       plt.ylabel('$u_{exact} - u$', fontsize=40)
#       plt.xlabel('$x$', fontsize=40)
#       plt.xticks(fontsize=30)
#       plt.yticks(fontsize=30)
#       plt.title('Errors  (%d points,   T = 0.075)'%pts[0], fontsize=40)
#       plt.plot(a, b, '--', linewidth=4, label= convection_type[i])

        temp = zip(a, b)

        with open('smooth-degree%d-'% (degree) + convection_type[i] + '-%d.csv'% pts[0], "w") as output:
           writer = csv.writer(output, lineterminator='\n', delimiter=',')
           for val in temp:
               writer.writerow(val)

#    plt.plot(a, c, linewidth=3, label='Exact')
#    leg = plt.legend(bbox_to_anchor=(1.00, 0.0), fontsize=40, loc='lower left', fancybox=True)
#    plt.tight_layout(rect=[0, 0, 0.65, 0.65])
#    plt.savefig('aliasing-bspline-error-%ddegree-%d.eps'% (degree, pts[0]), transparent=True, dpi=300)
#    plt.show()


def collocation_burgers_turbulence(degree, pts, p_deriv, convective_type):
    from periodic_derivatives import periodic_init
    from exact_periodic_burgers import exact_burgers
    from smooth_spline import smooth_spline_matrix_init
    from smooth_spline import smooth_spline_rhs_init
    from scipy.optimize import minimize_scalar
    from scipy.optimize import minimize

    # Set up formatting for the movie files
    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='burgers'), bitrate=1800)

    # save original number of points
    og_points = pts
    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2
    # time
    t = 0.0
    # runge kutta stage
    stage = 6
    # time
    t = 0.0
    # final time
    T = 0.05
    # time-step
    dt = 1.0e-5
    # diffusion
    nu = 5.0 * 1.0e-4

    u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)

    B = periodic_fix(ders[:, 0, :], degree)
    convection = periodic_fix(ders[:, 1, :], degree)
    diffusion = nu * periodic_fix(ders[:, 2, :], degree)

    u_n = np.zeros((stage, len(u)))

    # Initial data
    u_n[0, :] = ic_burgers_turbulence(pts)

    rk = runge_kutta(stage)

    print('%s, %d degree, %d points'%(convective_type, degree, pts))

    lamda = 1e-6
    A, firstcolumn, lastcolumn, fullB, P = \
            smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, p_deriv, 'periodic')

    BT = fullB.T

    # movie
    plt.rc('text', usetex=True)
    fig = plt.figure(figsize=(9, 4))
    axes = fig.add_subplot(211)
    uline, = axes.plot(u, u_n[0, :], '-', linewidth=2)
    frames = [u_n[0, :].copy()]
    tt = [0]
    num_plots = 200
    nplt = np.floor((T / num_plots) / dt)

    N = int(round(T / dt))

    LU = sla.lu_factor(B)
    smooth_LU = sla.lu_factor(A)

#    plt.rc('text', usetex=True)
#    fig = plt.figure(1, figsize=(9, 4))
#    ax = fig.add_subplot(111)
#    uline, = ax.plot(u, u_n[0, :], '-', linewidth=2)

    penalty = False

    shape, _ = BT.shape
    I = np.eye(shape)

    Prow = P[degree, :]

    for i in range(shape):
        P[i, :] = np.roll(Prow, -degree + i)

    u_ave = np.zeros_like(u_n[0, :])
    lam_all = np.zeros(N)
    for time in range(N):

#        if time == N-1:
#            penalty = True
#            lamda=1e-5
#            #lamda = minimize(objective, ((1.0e-5),), args=(fullB, I, P, u_n[0, :]), \
#            #        method='BFGS').x
#            #if (lamda <0.0 or lamda>1.0):
#            #    lamda = 0.0
#
#            #lam_all[time] = lamda

        umax = 0.5 * sum(abs(u_n[0, :]))**2
        # runge kutta SSP
        for m in range(stage):
            t += (rk[m] * dt)

            if (np.isnan(u_n[m, :]).any() or np.isinf(u_n[m, :]).any()):
                KE_time = np.linspace(0.0, time * dt, len(KE))
                return KE_time, KE

            if penalty == False:
                alpha = sla.lu_solve(LU, u_n[m, 1:])
            else:
                # smoothing spline
                smooth_LU = sla.lu_factor((BT @ fullB + lamda * P)[1:, 1:])
                rhs, bc = smooth_spline_rhs_init(u_n[m, :], firstcolumn, lastcolumn, BT, 'periodic')
                alpha = sla.lu_solve(smooth_LU, rhs)
                u_n[m, 1:] = B @ alpha

            if (convective_type == 'Advective'):
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                                    + 0.5 * dt * (diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              - 0.5 * dt * u_n[m, 1:] * (convection @ alpha) \
                              + 0.5 * dt * (diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Divergence'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * (
                                    - 0.5 * (convection @ alpha2) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              - 0.5 * (convection @ alpha2) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]

            elif (convective_type == 'Skew-symmetric'):
                alpha2 = sla.lu_solve(LU, u_n[m, 1:]**2)
                if ((m+1) != stage):
                    u_n[m + 1, 1:] = u_n[m, 1:]  \
                                    + 0.5 * dt * ( \
                                    -0.25 * (convection @ alpha2) \
                                    -0.50 * (u_n[m, 1:] * (convection @ alpha)) \
                                    + diffusion @ alpha)

                    u_n[m, 0] = u_n[m, -1]

                else:
                    u_n[0, 1:] = np.dot(rk[:m-1], u_n[:m-1, 1:]) + \
                                rk[-1] * (u_n[m, 1:] \
                              + 0.5 * dt * (
                              -0.25 * (convection @ alpha2) \
                              -0.50 * (u_n[0, 1:] * (convection @ alpha)) \
                              + diffusion @ alpha))

                    u_n[0, 0] = u_n[0, -1]


        #umax2 = 0.5 * sum(abs(u_n[0, :]))**2
        #if umax2 > umax:
        #    penalty = True
        #    print('Time = %e, L_inf = %e' % (t, umax2))
        #else:
        #    penalty = False

        print('Time = %e, Norm = %e' % (t, la.norm(u_n[0, :], np.inf)))
        assert(la.norm(u_n[0, :], np.inf) < 1e2)

        # saving frames for video
        if np.mod(time, nplt) == 0:
            frames.append(u_n[0, :].copy())
            tt.append(t)

#        uline.set_data(u, u_n[0, :])
#        ax.set_ylabel('$u(x, t)$', fontsize=30)
#        ax.set_xlabel('$x$', fontsize=30)
#        plt.ylim(top=0.65)
#        plt.xticks(fontsize=20)
#        plt.yticks(fontsize=20)
#        ax.set_title('B-splines, points = 256', fontsize=40)

#        if (time == N-1):
#            plt.tight_layout()
#            plt.savefig('op-burgers.eps', transparent=True, dpi=300)

#            t = np.linspace(0.0, T, len(lam_all))
#            plt.figure(2, figsize=(9, 4))
#            plt.ylabel('$\\lambda$', fontsize=30)
#            plt.xlabel('$t$', fontsize=30)
#            plt.xticks(fontsize=20)
#            plt.yticks(fontsize=20)
#            plt.plot(t, lam_all, '-', linewidth=2)
#            plt.title('Value of $\\lambda$ vs $t$', fontsize=40)
#            plt.tight_layout()
#            plt.savefig('op-lamda-burgers.eps', transparent=True, dpi=300)
#            plt.show()
#            quit()
        u_ave = u_ave + u_n[0, :]
    u_ave = u_ave / N
    data = np.array([u, u_n[0, :]]).T
    with open("burgers.out", 'w+') as file:
            np.savetxt(file, data, fmt=['%20.10f','%19.10f'])

    # movie
    im_ani = matplotlib.animation.FuncAnimation(fig, plot_frame, frames=len(frames),
            fargs=(uline, axes, u, frames, tt, pts), interval=20)
    im_ani.save('burgers.mp4', writer=writer)

    return u, u_n[0, :]

def main():
    # discrete_KE(5)
    # collocation_kdv(5, 128)
    # collocation_burgers(9, 32, 2, 'Advective')
    collocation_burgers_compact(9, 32, 2, 'Advective')
    # collocation_burgers_turbulence(7, 512, 2, 'Skew-symmetric')
    # collocation_burgers_smooth_vector(9, 32, 2, 'Advective')


if __name__ == "__main__":
    main()
