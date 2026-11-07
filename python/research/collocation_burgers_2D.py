import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import matplotlib.pyplot as plt
from collocation import nonperiodic_init
from smooth_spline import *
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
import csv


def objective(lamda, B, I, P, u):
#    Px = np.kron(lamda[0] * I, P)
#    Py = np.kron(P, lamda[1] * I)
#    Pxy = Px + Py
#    uhat = B @ la.solve(B.T @ B + Pxy, B.T @ u)
#    V = la.norm(u - uhat, 2)**2
#    V = V / np.trace(np.eye(u.shape[0]) - B @ la.inv(B.T @ B + Pxy) @ B.T)**2

    A = np.kron(B.T @ B + lamda[0]*P, B.T @ B + lamda[1]*P)
    fullB = np.kron(B.T @ B, B.T @ B)
    uhat = fullB @ la.solve(A, fullB.T @ u)
    V = la.norm(u - uhat, 2)**2

    V = V / np.trace(np.eye(u.shape[0]) - (fullB @ la.inv(A) @ fullB.T)/u.shape[0])**2
    return V / u.shape[0]


def runge_kutta(stage):
    if (stage == 6):
        return np.array([1.0/9.0, 2.0/5.0, 0.0, 4.0/9.0, 0.0, 2.0/45.0])
    elif (stage == 8):
        return np.array([2.0/15.0, 2.0/7.0, 2.0/9.0, 4.0/15.0, 0.0,
                        4.0/45.0, 0.0, 1.0/315.0])
    elif (stage == 10):
        return np.array([71./525, 22./81, 4./15, 4./21, 2./27, 4./75,
                        0.0, 8./945, 0.0, 2./14175])


def f(string, nu, x, y, t):
    import math as ma

    temp = np.exp((- t - 4.0 * x + 4.0 * y) / (32.0 * nu))

    if (string == 'u'):
        return 0.75 - 0.25 / (1.0 + temp)
    elif (string == 'v'):
        return 0.75 + 0.25 / (1.0 + temp)
    elif (string == 'lapu'):
        return - temp**2 / (64.0 * nu**2 * (1.0 + temp)**3) + temp / (128 * nu**2 * (1.0 + temp)**2)
    elif (string == 'lapv'):
        return temp**2 / (64.0 * nu**2 * (1.0 + temp)**3) - temp / (128 * nu**2 * (1.0 + temp)**2)
    elif (string == 'scalar-heat'):
        return 0.5 * (ma.erf((x + y + 1.0) / (2.0 * np.sqrt(2.0 * nu * t))) - \
                      ma.erf((x + y - 1.0) / (2.0 * np.sqrt(2.0 * nu * t))))


def slice_matrix(matrix, row, column):
    index = np.ix_(row, column)
    return matrix[index]


def xboundary(string, nu, t, fullu, xbc_index):

    n = len(fullu)

    x, y = np.unravel_index(xbc_index, (n, n))
    return f(string, nu, fullu[x], fullu[y], t).flatten('F')


def yboundary(string, nu, t, fullu, ybc_index):

    n = len(fullu)

    x, y = np.unravel_index(ybc_index, (n, n))
    return f(string, nu, fullu[x], fullu[y], t).flatten('F')


def construct_full_solution(alpha, xbc, ybc, fullu):

    n = len(fullu)

    b = np.append(xbc[0:n], alpha)
    b = np.append(b, xbc[n:])

    k = 0
    for i in range(n-2):
        b = np.insert(b, (i+1)*n    , ybc[k])
        b = np.insert(b, (i+2)*n-1, ybc[k+1])
        k += 2
    return b


def construct_reduced_solution(u, all_bc_index):
    return np.delete(u, all_bc_index)


def reduce_matrix(matrix):

    A = matrix.copy()
    # delete first and last row of matrix
    A = np.delete(A, 0, 0)
    A = np.delete(A, -1, 0)

    # delete first and last column of matrix
    A = np.delete(A, 0, 1)
    A = np.delete(A, -1, 1)
    return A


def indexing(fullu, string):

    n = len(fullu)

    all_index = np.arange(n**2, dtype=int)

    xbc_index = np.zeros(2*n, dtype=int)
    ybc_index = np.zeros(2*(n-2), dtype=int)

    # xbc row/column index
    xbc_index[:n] = all_index[:n]
    xbc_index[-n:] = all_index[-n:]

    # ybc row/column index
    k = 0
    for i in range(1, n-1):
        ybc_index[k] = i * n
        ybc_index[k+1] = (i + 1) * n - 1
        k += 2

    if string == 'x-nonbc':
        col_index = np.setdiff1d(xbc_index, ybc_index)

        row_index = np.setdiff1d(all_index, xbc_index)
        row_index = np.setdiff1d(row_index, col_index)
        row_index = np.setdiff1d(row_index, ybc_index)

    elif string == 'y-nonbc':
        col_index = np.setdiff1d(ybc_index, xbc_index)

        row_index = np.setdiff1d(all_index, xbc_index)
        row_index = np.setdiff1d(row_index, col_index)
        row_index = np.setdiff1d(row_index, ybc_index)

    elif string == 'bc':
        row_index = xbc_index
        col_index = ybc_index

    elif string == 'rhs-operator':
        row_index = np.setdiff1d(all_index, xbc_index)
        row_index = np.setdiff1d(row_index, ybc_index)

        col_index = row_index

    return all_index, row_index, col_index


def collocation_burgers_2D(degree, pts, penalty=False, output=False):
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # dimensions
    dim = 2
    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2
    # viscosity
    nu = 1.0/80.0
    # runge kutta stage
    stage = 6
    # vector component (i.e. u, v or w)
    component = ['u', 'v']
    # final time
    T = 1.0
    # time-step
    dt = 1.0e-2
    # time
    t = 0.0

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullD_1D = ders[:, 1, :]
    fullL_1D = ders[:, 2, :]

    n = len(fullu)

    # creating vectors
    Linf = np.zeros(dim)
    a = np.zeros(dim)
    #sol = np.zeros(n**2)
    xbc = np.zeros((dim, 2*n))
    ybc = np.zeros((dim, 2*(n-2)))

    # create full 2D mesh
    x, y = np.meshgrid(fullu, fullu)

    # create full operators
    m, m = fullL_1D.shape
    n = m**2
    fullI_1D = np.eye(m)
    fullI = np.eye(n)
    fullP = np.zeros((n, n))
    fullB = np.zeros((n, n))
    fullDx = np.zeros((n, n))
    fullDy = np.zeros((n, n))
    fullLap = np.zeros((n, n))

    _ ,_ ,_ ,_ , fullP_1D = \
    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    fullP_1D = fullP_1D * (maxk/pts)**2
    for i in range(m):
        for j in range(m):
            fullP[(i*m):((i+1)*m), (j*m):((j+1)*m)] = fullP_1D[i, j] * fullP_1D
            fullB[(i*m):((i+1)*m), (j*m):((j+1)*m)] = fullB_1D[i, j] * fullB_1D
            fullDx[(i*m):((i+1)*m), (j*m):((j+1)*m)] = fullD_1D[i, j] * fullB_1D
            fullDy[(i*m):((i+1)*m), (j*m):((j+1)*m)] = fullB_1D[i, j] * fullD_1D
            fullLap[(i*m):((i+1)*m), (j*m):((j+1)*m)] = fullL_1D[i, j] * fullB_1D + fullB_1D[i, j] * fullL_1D

    fullBinv = la.inv(fullB)
    uconvection = fullDx @ fullBinv
    vconvection = fullDy @ fullBinv
    diffusion = nu * (fullLap @ fullBinv)

    _, xbc_index, ybc_index = indexing(fullu, 'bc')
    _, xrow_index, xcol_index = indexing(fullu, 'x-nonbc')
    _, yrow_index, ycol_index = indexing(fullu, 'y-nonbc')
    _, row_index, col_index = indexing(fullu, 'rhs-operator')

    # flattened index of all bc points
    all_bc_index = np.append(xbc_index, ybc_index)

    # bc operators
    xbc_uconvection_matrix = slice_matrix(uconvection, xrow_index, xcol_index)
    ybc_uconvection_matrix = slice_matrix(uconvection, yrow_index, ycol_index)

    xbc_vconvection_matrix = slice_matrix(vconvection, xrow_index, xcol_index)
    ybc_vconvection_matrix = slice_matrix(vconvection, yrow_index, ycol_index)

    xbc_diffusion_matrix = slice_matrix(diffusion, xrow_index, xcol_index)
    ybc_diffusion_matrix = slice_matrix(diffusion, yrow_index, ycol_index)

    # operators
    I = np.eye(m-2)
    P = fullP_1D[1:-1, 1:-1]
    #P = slice_matrix(fullP, row_index, col_index)
    B = slice_matrix(fullB, row_index, col_index)
    Dx = slice_matrix(fullDx, row_index, col_index)
    Dy = slice_matrix(fullDy, row_index, col_index)
    Lap = slice_matrix(fullLap, row_index, col_index)

    LU = sla.lu_factor(B)

    #del fullB_1D
    del fullD_1D
    del fullL_1D
    #del fullB
    del fullDx
    del fullDy
    del fullLap
    del fullBinv
    del uconvection
    del vconvection
    del diffusion

    # time-step vectors
    n = len(fullu)

    fullu_n = np.zeros((dim, stage, n**2))
    u_n = np.zeros((dim, stage, (n-2)**2))
    fullalpha = np.zeros((dim, n**2))
    rk = runge_kutta(stage)
    lamda = np.zeros(dim)

    # initialize u
    for d in range(dim):
        fullu_n[d, 0, :] = f(component[d], nu, x, y, t).flatten('F')
        u_n[d, 0, :] = construct_reduced_solution(fullu_n[d, 0, :], all_bc_index)

    # number of time-steps
    N = int(T / dt)

#    if penalty == True:
#        fullmatrix = fullB.T @ fullB + 1e-9 * fullP
#        matrix = slice_matrix(fullmatrix, row_index, col_index)
#
#        LU = sla.lu_factor(matrix)
#
#        LUx = LU
#        LUy = LU

    for time in range(N):

        # SSP Runge-Kutta
        for m in range(stage):
            t += (rk[m] * dt)

            for d in range(dim):
                # x boundary
                xbc[d, :] = xboundary(component[d], nu, t, fullu, xbc_index)

                # y boundary
                ybc[d, :] = yboundary(component[d], nu, t, fullu, ybc_index)

                #if penalty == True:
                #    alpha = sla.lu_solve(LU, B.T @ u_n[d, 0, :])
                #else:
                #    alpha = sla.lu_solve(LU, u_n[d, 0, :])

                alpha = sla.lu_solve(LU, u_n[d, 0, :])

                #### FORM BCS
                bc = - u_n[0, m, :] * (xbc_uconvection_matrix @ xbc[d, :]) \
                     - u_n[1, m, :] * (xbc_vconvection_matrix @ xbc[d, :]) \
                     + (xbc_diffusion_matrix @ xbc[d, :])

                bc += - u_n[0, m, :] * (ybc_uconvection_matrix @ ybc[d, :]) \
                      - u_n[1, m, :] * (ybc_vconvection_matrix @ ybc[d, :]) \
                      + (ybc_diffusion_matrix @ ybc[d, :])

                bc = 0.5 * dt * bc

                # (m-stage, m-1 order) Explict SSP Runge-Kutta
                if ((m+1) != stage):
                    u_n[d, m + 1, :] = u_n[d, m, :]                             \
                                     - 0.5 * dt * (u_n[0, m, :] * (Dx @ alpha)) \
                                     - 0.5 * dt * (u_n[1, m, :] * (Dy @ alpha)) \
                                     + 0.5 * dt * (nu * Lap @ alpha)            \
                                     + bc
                else:
                    u_n[d, 0, :] = np.dot(rk[:m-1], u_n[d, :m-1, :])        \
                                 + rk[-1] * (u_n[d, m, :]                   \
                                 - 0.5 * dt * (u_n[0, m, :] * (Dx @ alpha)) \
                                 - 0.5 * dt * (u_n[1, m, :] * (Dy @ alpha)) \
                                 + 0.5 * dt * (nu * Lap @ alpha)            \
                                 + bc)

        if penalty == True:

            #fullu_n[0, 0, :] = construct_full_solution(u_n[0, 0, :], xbc[0, :], ybc[0, :], fullu)

            #res = minimize(objective, (1e-6, 1e-6), method='Nelder-Mead', options={'maxiter': 1000, \
            #       'xatol': 1e-8, 'fatol': 1e-8}, args=(fullB_1D, fullI_1D, fullP_1D, fullu_n[0, 0, :]))

            #if res.success == True and time != 0:
            #    if (time > 20 and res.x[0] > 1e-4 and res.x[1] > 1e-4):
            #        lamda1 = 1.6e-5
            #        lamda2 = 1.6e-5
            #    elif (res.x[0] < 0.0 or res.x[1] < 0.0):
            #        lamda1 = 0.0
            #        lamda2 = 0.0
            #    else:
            #        lamda1 = res.x[0]
            #        lamda2 = res.x[1]
            #else:
            #    lamda1 = 0.0
            #    lamda2 = 0.0

            #print(res.success, res.message, lamda1, lamda2)

            ## vector
            if time < 20:
                lamda1 = 0.0; lamda2 = 0.0
            else:
                lamda1 = 1.5e-5; lamda2 = 1.5e-5

            ### scalar
            #if time < 20:
            #    lamda1 = 0.0; lamda2 = 0.0
            #elif 20 <=time <= 60:
            #    lamda1 = 1.5e-5; lamda2 = 1.5e-5
            #else:
            #    lamda1 = 1e-4; lamda2 = 1.e-4

            A = np.kron(fullB_1D.T @ fullB_1D + lamda1*fullP_1D, fullB_1D.T @ fullB_1D + lamda2*fullP_1D)

            for d in range(dim):

                fullu_n[d, 0, :] = construct_full_solution(u_n[d, 0, :], xbc[d, :], ybc[d, :], fullu)
                fullu_n[d, 0, :] = fullB @ la.solve(A, fullB.T @ fullu_n[d, 0, :])
                u_n[d, 0, :] = construct_reduced_solution(fullu_n[d, 0, :], all_bc_index)

        # solution
        for d in range(dim):
            # numerical solution
            fullu_n[d, 0, :] = construct_full_solution(u_n[d, 0, :], xbc[d, :], ybc[d, :], fullu)

            # exact solution
            sol = f(component[d], nu, x, y, t).flatten('F')

            # L_inf error
            Linf[d] = la.norm(sol - fullu_n[d, 0, :], np.inf)
            assert(Linf[d] < 1e2)

            # location of max error
            a[d] = np.argmax(abs(sol - fullu_n[d, 0, :]))
            if np.intersect1d(int(a[d]), all_bc_index).size != 0:
                string = 'boundary'
            else:
                string = 'interior'

        print('Time = %f, uLinf = %e, vLinf = %e,  Points = %d, max error = %s' % (t, Linf[0], Linf[1], pts, string))

    #res = minimize(objective, (1e-5, 1e-5), method='BFGS', args=(fullB, fullI_1D, fullP_1D, fullu_n[d, 0, :]))
    #print(res.success)
    #print(res.x)
    #lamda1 = res.x[0]
    #lamda2 = res.x[1]

    #Px = np.kron(lamda1 * fullI_1D, fullP_1D)
    #Py = np.kron(fullP_1D, lamda2 * fullI_1D)
    #Pxy = Px + Py
    #fullu_n[0, 0, :] = fullB @ la.solve(fullB.T @ fullB + Pxy, fullB.T @ fullu_n[0, 0, :])

    #print(la.norm(sol - fullu_n[0, 0, :], np.inf))

    x, y = np.meshgrid(fullu, fullu)
    zx = np.reshape(fullu_n[0, 0, :], (pts, pts), order='F')
    zy = np.reshape(fullu_n[1, 0, :], (pts, pts), order='F')
    ax.plot_surface(x, y, zx)
    plt.show()

    if output==True and penalty==True:
        with open('velocity-smooth-collocation-burgers-2D-degree%d.csv'%degree, "w") as output:
            writer = csv.writer(output, lineterminator='\n', delimiter=',')
            for j in range(len(x)):
                for i in range(len(x)):
                    writer.writerow([x[i,j], y[i,j], zx[i,j], zy[i,j]])

    elif output==True and penalty==False:
        with open('velocity-original-collocation-burgers-2D-degree%d.csv'%degree, "w") as output:
            writer = csv.writer(output, lineterminator='\n', delimiter=',')
            for j in range(len(x)):
                for i in range(len(x)):
                    writer.writerow([x[i,j], y[i,j], zx[i,j], zy[i,j]])

    return Linf[0], Linf[1]


def test_collocation_burgers_2D(degree):

    uerror = []
    verror = []
    pts = np.array([30, 40])[::-1]

    for i in range(len(pts)):
        a, b = collocation_burgers_2D(degree, pts[i])

        uerror.append(a)
        verror.append(b)

        temp = zip(pts, uerror, verror)

        with open('collocation-burgers-2D-degree%d.csv'% degree, "w") as output:
            writer = csv.writer(output, lineterminator='\n', delimiter=',')
            for val in temp:
                writer.writerow(val)


def main():
    collocation_burgers_2D(7, 20, penalty=False, output=True)
    # test_collocation_burgers_2D(5)


if __name__ == "__main__":
    main()
