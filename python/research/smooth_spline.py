import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from wavenumber_periodic import periodic_fix
from periodic_derivatives import periodic_init
from collocation import nonperiodic_init


def noisy_sin(u):
    y = np.zeros_like(u)
    for i in range(len(u)):
        if 0.25 *  np.max(u) < u[i] < 0.5 * np.max(u):
            y[i] = np.sin(8 * u[i]) + (2.0 * (np.random.random(1) - 0.5)) * 0.3
        else:
            y[i] = np.sin(8 * u[i])
    return y


def noisy_exp(u):
    y = np.exp(-10.0 * (u - np.mean(u))**2)
    for i in range(len(y)):
        if abs(y[i]) > 1e-2:
            y[i] += np.random.random(1) * 0.3
    return y


def noisy_step(u):
    y = np.zeros_like(u)
    for i in range(len(u)):
        if u[i] < np.mean(u):
            y[i] = 2.0
        else:
            y[i] = 1.0
        #if 0.25 *  np.max(u) < u[i] < 0.5 * np.max(u):
        #    y[i] = y[i] + (2.0 * (np.random.random(1) - 0.5)) * 0.3
    return y


def square_wave(u):
    return np.sign(np.sin(2.0 * u))


def sawtooth(u):
    y = np.zeros_like(u)
    for i in range(len(u)):
        if 0.0 <= u[i] <= np.pi:
            y[i] = 3.0 * u[i]
        else:
            y[i] = 3.0 * (u[i] - 2.0*np.pi)
    return y


def noisy_tanh(u):
    y = np.zeros_like(u)
    for i in range(len(u)):
        if u[i] < np.mean(u):
            y[i] = np.tanh(500.0 * (u[i] - np.mean(u)))
        else:
            y[i] = np.tanh(500.0 * (u[i] - np.mean(u)))
    return y


def smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, deriv, bc_type, matrix_type='reduced'):
    from scipy.integrate import quad
    from galerkin import B_prime
    from knot_multiplicity import knot_multiplicity
    from knot_multiplicity import periodic_knot_multiplicity

    if bc_type == 'periodic':
        pts = pts - (degree - 1)
        u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)
        #_, _, U2 = nonperiodic_init(mink, maxk, pts, degree, deriv)

        U = np.linspace(mink, maxk, pts)
        U2 = knot_multiplicity(degree, U)
        U = periodic_knot_multiplicity(degree, U)

    elif bc_type == 'nonperiodic':
        _, ders, U = nonperiodic_init(mink, maxk, pts, degree, deriv)
        U2 = U.copy()
    else:
        print('Invalid bc_type. Choose \'periodic\' or \'non-periodic\'')
        quit()

    size, n_splines = ders[:, 0, :].shape
    A = np.zeros((n_splines, n_splines))
    D = np.zeros((n_splines, n_splines))
    x = np.zeros((deriv + 1, degree + 1))

    # position of B_i in respective subintervals
    a = np.arange(degree + 1)[::-1]

    xi = np.zeros(degree + 2)

    row, col = D.shape
    # quad rule with B-splines
    for i in range(row):

        # coordinates of subintervals
        for j in range(len(xi)):
            xi[j] = U2[i + j]

        for band in range(-degree, degree + 1):
            for j in range(len(a) - abs(band)):

                if (band >= 0 and (i + band) < col):
                    k = j + band

                    # (B_j'', B_i'')
                    D[i, i + band] += quad(B_prime, xi[k], xi[k+1], args=(a[k], x, deriv, n_splines, U, degree, band))[0]

                elif (band < 0 and (i + band) >= 0):
                    # (B_j'', B_i'')
                    D[i, i + band] += quad(B_prime, xi[j], xi[j+1], args=(a[j], x, deriv, n_splines, U, degree, band))[0]

    row, col = A.shape

    if bc_type == 'periodic':
        I = np.zeros((row, col))
        periodic_row = np.zeros(row)
        periodic_row[:degree] = B[0, :degree]

        for i in range(row):
            I[i, :] = np.roll(periodic_row, i)

    elif bc_type == 'nonperiodic':
        I = B

    A = (I.T @ I) + (lamda * D)

    if bc_type == 'nonperiodic':
        if matrix_type == 'reduced':
            # delete first/last row
            A = np.delete(A, 0, 0)
            A = np.delete(A, -1, 0)

            # save first/last column
            firstcolumn = A[:, 0]
            lastcolumn = A[:, -1]

            # delete first/last column
            A = np.delete(A, 0, 1)
            A = np.delete(A, -1, 1)
        else:
            firstcolumn = 0.0
            lastcolumn = 0.0

    elif bc_type == 'periodic':
        # delete first row
        A = np.delete(A, 0, 0)

        # save first/last column
        firstcolumn = A[:, 0]
        lastcolumn = A[:, -1]

        # delete first row
        A = np.delete(A, 0, 1)

    return A, firstcolumn, lastcolumn, I, D


def smooth_spline_rhs_init(y, firstcolumn, lastcolumn, A, bc_type):

    if bc_type == 'nonperiodic':

        # rhs and bc
        bc0 = y[0]
        bc1 = y[-1]
        b = y[1:-1].copy()

        # place bc to rhs
        b = b - bc0 * firstcolumn
        b = b - bc1 * lastcolumn

        # delete first/last row
        A = np.delete(A, 0, 0)
        A = np.delete(A, -1, 0)

        # delete first/last column
        A = np.delete(A, 0, 1)
        A = np.delete(A, -1, 1)

        b = A @ b

        return b, np.array([bc0, bc1])
    elif bc_type == 'periodic':

        # rhs and bc
        bc0 = y[0]
        b = y[1:].copy()

        # place bc to rhs
        b = b - bc0 * firstcolumn

        # delete first row/column
        A = np.delete(A, 0, 0)
        A = np.delete(A, 0, 1)

        b = A @ b

        return b, np.array([bc0])

def test_smooth_burgers():
    from collocation_kdv import objective
    from scipy.optimize import minimize

    mink = 0.0
    maxk = 2.0*np.pi
    lamda =1.5e-6
    degree = 5

    u, un = np.loadtxt("results/burgers1d/gibbs/bspline-burgers.out", skiprows=0, unpack=True)
    pts = len(u)
    u, ders, _ = nonperiodic_init(mink, maxk, pts, degree, 2)
    B = ders[:, 0, :]
    I = np.eye(B.shape[0])

    _, _, _,_ , P = smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, 2, 'nonperiodic')

    res = minimize(objective, (1e-5,), method='Nelder-Mead', options={'maxiter': 5000, \
           'xatol': 1e-12, 'fatol': 1e-12}, args=(B, I, P, un))
    lamda = res.x
    print(res.success, res.message, res.x)

    m, n = B.shape

    un = la.solve(B.T @ B + lamda * P, B.T @ un)
    un = B @ un

    data = np.array([u, un]).T
    with open("results/burgers1d/gibbs/filter-burgers.out", 'w+') as datafile_id:
            np.savetxt(datafile_id, data, fmt=['%20.10f','%19.10f'])

    plt.rc('text', usetex=True)
    fig = plt.figure(1, figsize=(9, 4))
    plt.plot(u, un, linewidth=2)
    plt.ylabel('$u(x, t)$', fontsize=30)
    plt.xlabel('$x$', fontsize=30)
    plt.ylim(top=0.65)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('B-spline, points = 256', fontsize=40)
    plt.tight_layout()
    plt.show()



def test_smooth_spline(degree, pts, bc_type, test):
    from collocation_kdv import objective
    from scipy.optimize import minimize

    mink = 0.0
    maxk = 2.0*np.pi

    plt.rc('text', usetex=True)
    fig1 = plt.figure(1, figsize=(18, 14))
    #fig2 = plt.figure(2, figsize=(18, 14))
    ax1 = fig1.add_subplot(111)
    #ax2 = fig2.add_subplot(111)
    #ax2.set_title('B-spline coefficients (noisy $-$ filtered)', fontsize=70)

    if test == 'exp':
        all_lamda = [0.002]
        #all_lamda = [0.001, 0.002, 0.003, 0.004, 0.005]
        function = noisy_exp

        ax1.set_title('$\exp(-10 (x-\pi)^2)$ + noise', fontsize=70)
    elif test == 'sin':
        #lamda = 0.02
        all_lamda = [1.0e-8]
        function = noisy_sin

        ax1.set_title('$\sin(8x)$ + noise', fontsize=70)
        #ax2.set_title('B-spline coefficients (noisy $-$ filtered)', fontsize=70)

    elif test == 'squarewave':
        all_lamda = [0.002, 0.0]
        function = square_wave

        ax1.set_title('Square Waves', fontsize=50)
        #ax2.set_title('B-spline coefficients (noisy $-$ filtered)', fontsize=70)

    elif test == 'sawtooth':
        all_lamda = [0.002, 0.0]
        function = sawtooth

        ax1.set_title('Sawtooth', fontsize=50)
        #ax2.set_title('B-spline coefficients (noisy $-$ filtered)', fontsize=70)

    elif test == 'step':
        all_lamda = [0.002]
        function = noisy_step

        ax1.set_title('Step function + noise', fontsize=70)
        #ax2.set_title('B-spline coefficients (noisy $-$ filtered)', fontsize=70)

    elif test == 'tanh':
        all_lamda = [1.0e-8]
        function = noisy_tanh

        ax1.set_title('$\\tanh(500(x-\pi))$ + noise', fontsize=70)
        #ax2.set_title('B-spline coefficients (noisy $-$ filtered)', fontsize=70)

    for lam in range(len(all_lamda)):
        lamda = all_lamda[lam]
        if bc_type == 'periodic':
            u, ders, _ = periodic_init(mink, maxk, pts, degree, 2)
            B = periodic_fix(ders[:, 0, :], degree)
            rhs = function(u)
        elif bc_type == 'nonperiodic':
            u, ders, _ = nonperiodic_init(mink, maxk, pts, degree, 2)
            B = ders[:, 0, :]
            B_og = B.copy()
            rhs = function(u)

        if test == "squarewave":
            rhs[0] = 1.0

        _, _, _,_ , P = smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, 2, bc_type)

        m, n = B.shape
        I = np.eye(m)

        if lamda != 0:
            res = minimize(objective, (1e-6,), method='Nelder-Mead', options={'maxiter': 1000, \
                   'xatol': 1e-10, 'fatol': 1e-10}, args=(B, I, P, rhs))
            lamda1 = res.x
            print(res.success, res.message, res.x)

            if test == 'squarewave':
                lamda1 = lamda1 * 1e-6
            elif test == 'sawtooth':
                lamda1 = lamda1 * 5.6832942341 * 1e-6
        lamda1 = 1e-4
        alpha_filtered = la.solve(B.T @ B + lamda1 * P, B.T @ rhs)

        if bc_type == 'periodic':

            yfiltered = np.insert(yfiltered, 0, bc[0])

            yfiltered[1:] = B @ yfiltered[1:]
            yfiltered[-1] = yfiltered[0]

        elif bc_type == 'nonperiodic':
            alpha_og = la.solve(B, rhs)
            yfiltered = B @ alpha_filtered
            alpha_filtered[:3*degree] = alpha_og[:3*degree]
            alpha_filtered[-3*degree:] = alpha_og[-3*degree:]

        if lam == 0:
            ax1.plot(u, rhs, 'o', ms=9, label='data')
            ax1.plot(u, rhs, '--', linewidth=5, label='B-spline')
        if lamda == 0:
            ax1.plot(u, rhs, '--', linewidth=5, label='B-splines, $B^5$')
            #ax2.plot(u, alpha_og, '--', linewidth=6, label='noisy')

            if test == 'sin':
                exact = np.sin(8.0 * u)
                error_bspline = np.max(np.abs(exact - rhs))
                error_filtered = np.max(np.abs(exact - yfiltered))
                print("Error from B-spline %e"%error_bspline)
                print("Error from filter %e"%error_filtered)
            elif test == 'exp':
                exact = np.exp(-10.0 * (u - np.mean(u))**2)
                error_bspline = np.max(np.abs(exact - rhs))
                error_filtered = np.max(np.abs(exact - yfiltered))
                print("Error from B-spline %e"%error_bspline)
                print("Error from filter %e"%error_filtered)
            elif test == 'squarewave':
                y1 = (4/np.pi)*np.sin(2.0 * u)
                y3 = y1 + (4/(3*np.pi))*np.sin(2.0 * 3*u)
                y5 = y3 + (4/(5*np.pi))*np.sin(2.0 * 5*u)
                ax1.plot(u, y5, '--', linewidth=5, label='Fourier Series, $S_5$')
            elif test == 'sawtooth':
                y1 = np.sin(np.pi/np.pi*u)
                y2 = y1 - 1/2.0*np.sin(2.0*np.pi/np.pi*u)
                y3 = y2 + 1/3.0*np.sin(3.0*np.pi/np.pi*u)
                y4 = y3 - 1/4.0*np.sin(4.0*np.pi/np.pi*u)
                y5 = y4 + 1/5.0*np.sin(5.0*np.pi/np.pi*u)
                ax1.plot(u, 6*np.pi/np.pi*y5, '--', linewidth=5, label='Fourier Series, $S_5$')

        else:
            #ax1.plot(u, yfiltered, linewidth=5, label='smooth, $\lambda^2$=%.2e'%abs(lamda1))
            ax1.plot(u, yfiltered, linewidth=5, label='smooth')
            #ax2.plot(u, alpha_og - alpha_filtered, '--', linewidth=6, label='$\lambda$=%.2e'%lamda)

            ax1.tick_params(axis='both', labelsize=40)
            #ax1.legend(fontsize=40, loc='upper right')
            #ax2.tick_params(axis='both', labelsize=40)
            #ax2.legend(fontsize=40, loc='upper right')

    #leg = plt.legend(bbox_to_anchor= (1.02, 0), loc='lower left', fontsize=35)
    #leg.get_frame().set_alpha(0.3)
    #plt.tight_layout(rect=[0, 0, 0.55, 0.65])

    plt.legend(loc='upper right', fontsize=35)
    plt.tight_layout()

    plt.xlim(mink, maxk)
    fig1.savefig('u-%s-noise.eps'%test, transparent=True, dpi=300)
    #fig2.savefig('alpha-%s-noise.eps'%test, transparent=True, dpi=300)
    plt.show()


def main():
    test_smooth_spline(5, 100, 'nonperiodic', 'sin')
    #test_smooth_burgers()


if __name__ == "__main__":
    main()
