import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import matplotlib.pyplot as plt
from collocation import nonperiodic_init
from periodic_derivatives import periodic_init
from wavenumber_periodic import periodic_fix
from collocation_ns_common import *


def linear_system(xLU, yLU, rho, u, v, T, tau,
                  alpha, alpha_T, alpha_tau,
                  alpha_rho_u, alpha_rho_v, alpha_rho_T,
                  alpha_conv_T, alpha_conv_div):

    # linear solve
    x, _ = xLU.shape
    y, _ = yLU.shape

    alpha_T[:]              = (yLU.solve(xLU.solve(np.reshape(T, (x, y), 'C')).T).T).flatten('C')
    alpha[0, :]             = (yLU.solve(xLU.solve(np.reshape(u, (x, y), 'C')).T).T).flatten('C')
    alpha[1, :]             = (yLU.solve(xLU.solve(np.reshape(v, (x, y), 'C')).T).T).flatten('C')
    alpha_rho_u[:]          = (yLU.solve(xLU.solve(np.reshape(rho * u, (x, y), 'C')).T).T).flatten('C')
    alpha_rho_v[:]          = (yLU.solve(xLU.solve(np.reshape(rho * v, (x, y), 'C')).T).T).flatten('C')
    alpha_rho_T[:]          = (yLU.solve(xLU.solve(np.reshape(rho * T, (x, y), 'C')).T).T).flatten('C')
    alpha_conv_div[0, 0, :] = (yLU.solve(xLU.solve(np.reshape(rho * u**2, (x, y), 'C')).T).T).flatten('C')
    alpha_conv_div[1, 1, :] = (yLU.solve(xLU.solve(np.reshape(rho * v**2, (x, y), 'C')).T).T).flatten('C')
    alpha_conv_div[0, 1, :] = (yLU.solve(xLU.solve(np.reshape(rho * u * v, (x, y), 'C')).T).T).flatten('C')
    alpha_conv_div[1, 0, :] = alpha_conv_div[0, 1, :]
    alpha_conv_T[0, :]      = (yLU.solve(xLU.solve(np.reshape(rho * T * u, (x, y), 'C')).T).T).flatten('C')
    alpha_conv_T[1, :]      = (yLU.solve(xLU.solve(np.reshape(rho * T * v, (x, y), 'C')).T).T).flatten('C')
    alpha_tau[0, 0, :]      = (yLU.solve(xLU.solve(np.reshape(tau[0, 0, :], (x, y), 'C')).T).T).flatten('C')
    alpha_tau[0, 1, :]      = (yLU.solve(xLU.solve(np.reshape(tau[0, 1, :], (x, y), 'C')).T).T).flatten('C')
    alpha_tau[1, 0, :]      = alpha_tau[0, 1, :]
    alpha_tau[1, 1, :]      = (yLU.solve(xLU.solve(np.reshape(tau[1, 1, :], (x, y), 'C')).T).T).flatten('C')
    return


def periodic_boundary_condition(u_n_ybc, T_ybc, rho_T_ybc, rho_u_ybc, rho_v_ybc,
                                conv_div_ybc, rho_T_u_ybc, rho_T_v_ybc,
                                rkt, fullu, ybc_index, function, primitive, test, rho_ybc=None):
    if test == 'mms':
        # y boundary

        # velocity u
        u_n_ybc[0, :] = yboundary(function[1], rkt, fullu, ybc_index)
        # velocity v
        u_n_ybc[1, :] = yboundary(function[2], rkt, fullu, ybc_index)
        # temperature T
        T_ybc[:] = yboundary(function[3], rkt, fullu, ybc_index)
        # rho * T
        rho_T_ybc[:] = yboundary(primitive[5], rkt, fullu, ybc_index)
        # rho * u
        rho_u_ybc[:] = yboundary(primitive[0], rkt, fullu, ybc_index)
        # rho * v
        rho_v_ybc[:] = yboundary(primitive[1], rkt, fullu, ybc_index)
        # rho * u * u
        conv_div_ybc[0, 0, :] = yboundary(primitive[2], rkt, fullu, ybc_index)
        # rho * u * v
        conv_div_ybc[0, 1, :] = yboundary(primitive[3], rkt, fullu, ybc_index)
        # rho * v * u
        conv_div_ybc[1, 0, :] = conv_div_ybc[0, 1, :]
        # rho * v * v
        conv_div_ybc[1, 1, :] = yboundary(primitive[4], rkt, fullu, ybc_index)
        # rho * T * u
        rho_T_u_ybc[:] = yboundary(primitive[6], rkt, fullu, ybc_index)
        # rho * T * v
        rho_T_v_ybc[:] = yboundary(primitive[7], rkt, fullu, ybc_index)

    elif test == 'channel':

        # velocity u
        u_n_ybc[0, :] = 0.0
        # velocity v
        u_n_ybc[1, :] = 0.0
        # temperature T
        T_ybc[:] = 1.0
        # rho * T
        rho_T_ybc[:] = rho_ybc * T_ybc
        # rho * u
        rho_u_ybc[:] = rho_ybc * u_n_ybc[0, :]
        # rho * v
        rho_v_ybc[:] = rho_ybc * u_n_ybc[1, :]
        # rho * u * u
        conv_div_ybc[0, 0, :] = rho_ybc * u_n_ybc[0, :]**2
        # rho * u * v
        conv_div_ybc[0, 1, :] = rho_ybc * u_n_ybc[0, :] * u_n_ybc[1, :]
        # rho * v * u
        conv_div_ybc[1, 0, :] = conv_div_ybc[0, 1, :]
        # rho * v * v
        conv_div_ybc[1, 1, :] = rho_ybc * u_n_ybc[1, :]**2
        # rho * T * u
        rho_T_u_ybc[:] = rho_ybc * T_ybc * u_n_ybc[0, :]
        # rho * T * v
        rho_T_v_ybc[:] = rho_ybc * T_ybc * u_n_ybc[1, :]
    return


def continuity_bc(cont_bc, ybc_DBinv, rho_u_ybc, rho_v_ybc):

    # d(rhou)dx + d(rhov)dy
    cont_bc[:] = ybc_DBinv[0] @ rho_u_ybc + ybc_DBinv[1] @ rho_v_ybc
    return


def momentum_bc(conv_bc, Lap_bc, rho_T_bc, tau, ybc_tau,
                fullu, xLU, yfullLU, tempfullU, xboundary,
                tempU, tempRho, u_n_ybc, rho_T_ybc, conv_div_ybc, rho_u_ybc, rho_v_ybc,
                all_bc_index, ybc_index, fullDx, fullDy, alpha_full_u, mu,
                ybc_DBinv):

    x_periodic_velocity_gradient_bc(tau, ybc_tau,
        tempU, u_n_ybc, fullu, xboundary, xLU, yfullLU, tempfullU,
        all_bc_index, ybc_index, fullDx, fullDy, alpha_full_u, mu)

    for d in range(2):
        # rho * u * dudx + rho * v * dudy
        conv_bc[d, :] = tempRho * (tempU[0, :] * (ybc_DBinv[0] @ u_n_ybc[d, :]) + \
                                   tempU[1, :] * (ybc_DBinv[1] @ u_n_ybc[d, :]))

        # d(rho * u * u ) / dx + d(rho * u * v) / dy
        conv_bc[d, :] += ybc_DBinv[0] @ conv_div_ybc[d, 0, :] + ybc_DBinv[1] @ conv_div_ybc[d, 1, :]

        # -u * (d(rho * u) / dx + d(rho * v) / dy)
        conv_bc[d, :] += - tempU[d, :] * (ybc_DBinv[0] @ rho_u_ybc + ybc_DBinv[1] @ rho_v_ybc)

        # dti1/dx + dti2/dy
        Lap_bc[d, :] = ybc_DBinv[0] @ ybc_tau[d, 0, :] + ybc_DBinv[1] @ ybc_tau[d, 1, :]

        # d(rho * T)/dx_i
        rho_T_bc[d, :] = ybc_DBinv[d] @ rho_T_ybc
    return


def energy_bc(Tconv_bc, TLap_bc, Tdudx_bc, Tdivu_bc,
              tempU, tempT, tempRho,
              u_n_ybc, T_ybc, rho_u_ybc, rho_v_ybc, rho_T_u_ybc, rho_T_v_ybc,
              ybc_DBinv, ybc_D2Binv):

    # rho * u * dTdx + rho * v * dTdy
    Tconv_bc[:] = tempRho * (tempU[0, :] * (ybc_DBinv[0] @ T_ybc) + \
                             tempU[1, :] * (ybc_DBinv[1] @ T_ybc))

    # d(rho * T * u ) / dx + d(rho * T * v) / dy
    Tconv_bc[:] += ybc_DBinv[0] @ rho_T_u_ybc + ybc_DBinv[1] @ rho_T_v_ybc

    # -T * (d(rho * u) / dx + d(rho * v) / dy)
    Tconv_bc[:] += - tempT * (ybc_DBinv[0] @ rho_u_ybc + ybc_DBinv[1] @ rho_v_ybc)

    # Dxx T + Dyy T
    TLap_bc[:] = ybc_D2Binv[0] @ T_ybc + ybc_D2Binv[1] @ T_ybc

    # dudx bc
    Tdudx_bc[0, 0, :]= ybc_DBinv[0] @ u_n_ybc[0, :]
    # dudy bc
    Tdudx_bc[0, 1, :]= ybc_DBinv[1] @ u_n_ybc[0, :]
    # dvdx bc
    Tdudx_bc[1, 0, :]= ybc_DBinv[0] @ u_n_ybc[1, :]
    # dvdy bc
    Tdudx_bc[1, 1, :]= ybc_DBinv[1] @ u_n_ybc[1, :]

    # T * (dudx + dvdy)
    Tdivu_bc[:] = tempT * (Tdudx_bc[0, 0, :] + Tdudx_bc[1, 1, :])
    return


def boundary_condition(rkt, fullu, xbc_index, ybc_index,
                       function, primitive, source,
                       u_n_xbc, T_xbc, rho_T_xbc, rho_u_xbc, rho_v_xbc,
                       conv_div_xbc, rho_T_u_xbc, rho_T_v_xbc,
                       u_n_ybc, T_ybc, rho_T_ybc, rho_u_ybc, rho_v_ybc,
                       conv_div_ybc, rho_T_u_ybc, rho_T_v_ybc):

    # x boundary

    # velocity u
    u_n_xbc[0, :] = xboundary(function[1], rkt, fullu, xbc_index)
    # velocity v
    u_n_xbc[1, :] = xboundary(function[2], rkt, fullu, xbc_index)
    # temperature T
    T_xbc[:] = xboundary(function[3], rkt, fullu, xbc_index)
    # rho * T
    rho_T_xbc[:] = xboundary(primitive[5], rkt, fullu, xbc_index)
    # rho * u
    rho_u_xbc[:] = xboundary(primitive[0], rkt, fullu, xbc_index)
    # rho * v
    rho_v_xbc[:] = xboundary(primitive[1], rkt, fullu, xbc_index)
    # rho * u * u
    conv_div_xbc[0, 0, :] = xboundary(primitive[2], rkt, fullu, xbc_index)
    # rho * u * v
    conv_div_xbc[0, 1, :] = xboundary(primitive[3], rkt, fullu, xbc_index)
    # rho * v * u
    conv_div_xbc[1, 0, :] = conv_div_xbc[0, 1, :]
    # rho * v * v
    conv_div_xbc[1, 1, :] = xboundary(primitive[4], rkt, fullu, xbc_index)
    # rho * T * u
    rho_T_u_xbc[:] = xboundary(primitive[6], rkt, fullu, xbc_index)
    # rho * T * v
    rho_T_v_xbc[:] = xboundary(primitive[7], rkt, fullu, xbc_index)

    # y boundary

    # velocity u
    u_n_ybc[0, :] = yboundary(function[1], rkt, fullu, ybc_index)
    # velocity v
    u_n_ybc[1, :] = yboundary(function[2], rkt, fullu, ybc_index)
    # temperature T
    T_ybc[:] = yboundary(function[3], rkt, fullu, ybc_index)
    # rho * T
    rho_T_ybc[:] = yboundary(primitive[5], rkt, fullu, ybc_index)
    # rho * u
    rho_u_ybc[:] = yboundary(primitive[0], rkt, fullu, ybc_index)
    # rho * v
    rho_v_ybc[:] = yboundary(primitive[1], rkt, fullu, ybc_index)
    # rho * u * u
    conv_div_ybc[0, 0, :] = yboundary(primitive[2], rkt, fullu, ybc_index)
    # rho * u * v
    conv_div_ybc[0, 1, :] = yboundary(primitive[3], rkt, fullu, ybc_index)
    # rho * v * u
    conv_div_ybc[1, 0, :] = conv_div_ybc[0, 1, :]
    # rho * v * v
    conv_div_ybc[1, 1, :] = yboundary(primitive[4], rkt, fullu, ybc_index)
    # rho * T * u
    rho_T_u_ybc[:] = yboundary(primitive[6], rkt, fullu, ybc_index)
    # rho * T * v
    rho_T_v_ybc[:] = yboundary(primitive[7], rkt, fullu, ybc_index)
    return


def convection_eval(DxBinv, DyBinv, rho, u, v, convection_type):

    if convection_type == 'convection divergence':
        # d(rho * phi)/dt + d(rho *phi * u_j)/dx_j
        return  np.array([DxBinv @ (rho * u * u) + DyBinv @ (rho * u * v),\
                          DxBinv @ (rho * v * u) + DyBinv @ (rho * v * v)])
    elif convection_type == 'convection advective':
        # rho * d(phi)/dt + rho * u_j * d(phi)/dx_j
        return  np.array([(rho * u) * (DxBinv @ u) + (rho * v) * (DyBinv @ u),\
                          (rho * u) * (DxBinv @ v) + (rho * v) * (DyBinv @ v)])
    elif convection_type == 'convection skew-symmetric':
        # rho d(phi)/dt + 1/2 * [rho * u_j * d(phi)/dx_j + d(rho * phi * u_j)/dx_j
        #                       - phi * d(rho * u_j)/dx_j]
        return  0.5 * np.array([((rho * u) * (DxBinv @ u) + (rho * v) * (DyBinv @ u) + \
                               DxBinv @ (rho * u * u) + DyBinv @ (rho * u * v) -       \
                               u * (DxBinv @ (rho * u) + DyBinv @ (rho * v))),         \
                               ((rho * u) * (DxBinv @ v) + (rho * v) * (DyBinv @ v) +  \
                               DxBinv @ (rho * v * u) + DyBinv @ (rho * v * v) -       \
                               v * (DxBinv @ (rho * u) + DyBinv @ (rho * v)))])
    else:
        print('No such convection_type: \'%s\''%convection_type)
        quit()


def test_convection_eval(convection, convection_type, x, y, L, t):

    # solution
    sol = np.array([f('u' + convection_type, x, y, L, t).flatten('F'), \
                    f('v' + convection_type, x, y, L, t).flatten('F')])

    # error
    Linf = np.zeros(len(convection))
    for i in range(len(convection)):
        Linf[i] = la.norm(convection[i] - sol[i], np.inf)

    # its just the way they index the numpy mesh
    ny, nx = x.shape

    print(nx, ny, Linf[0], Linf[1])
    return Linf[0], Linf[1]


def viscous_operator(DxBinv, DyBinv, Dx=None, Dy=None):

    m, n = DyBinv.shape
    assert(m == n)
    m1, n1 = DxBinv.shape
    assert(m1 == n1)
    assert(m1 == n)

    viscous_operator = np.zeros((2*m, 2*m))

    if Dx is None and Dy is None:
        viscous_operator[0:m, 0:m] = DxBinv @ (4.0/3.0 * DxBinv) + DyBinv @ DyBinv
        viscous_operator[0:m, m: ] = DxBinv @ (-2.0/3.0 * DyBinv) + DyBinv @ DxBinv
        viscous_operator[m: , 0:m] = viscous_operator[0:m, m:]
        viscous_operator[m: , m: ] = DxBinv @ DxBinv + DyBinv @ (4.0/3.0 * DyBinv)
    else:
        viscous_operator[0:m, 0:m] = DxBinv @ (4.0/3.0 * Dx) + DyBinv @ Dy
        viscous_operator[0:m, m: ] = DxBinv @ (-2.0/3.0 * Dy) + DyBinv @ Dx
        viscous_operator[m: , 0:m] = viscous_operator[0:m, m:]
        viscous_operator[m: , m: ] = DxBinv @ Dx + DyBinv @ (4.0/3.0 * Dy)
    return viscous_operator


def x_periodic_velocity_gradient_bc(tau, ybc_tau,
        tempU, u_n_ybc, fullu, xboundary, xLU, yfullLU, tempfullU,
        all_bc_index, ybc_index, fullDx, fullDy, alpha_full_u, mu):

    nx = len(fullu[0])
    ny = len(fullu[1])

    x, _ = xLU.shape
    y, _ = yfullLU.shape

    for d in range(2):
        tempfullU[d, :] = construct_full_solution(tempU[d, :], None, u_n_ybc[d, :], fullu, xboundary)[:(nx*ny)]
        alpha_full_u[d, :] = (yfullLU.solve(xLU.solve(np.reshape(tempfullU[d, :], (x, y), 'C')).T).T).flatten('C')

    t11 = mu * (4.0/3.0 * (fullDx @ alpha_full_u[0, :]) - 2.0/3.0 * (fullDy @ alpha_full_u[1, :]))
    t12 = mu * (fullDy @ alpha_full_u[0, :] + fullDx @ alpha_full_u[1, :])
    t21 = t12
    t22 = mu * (4.0/3.0 * (fullDy @ alpha_full_u[1, :]) - 2.0/3.0 * (fullDx @ alpha_full_u[0, :]))

    # tau bc
    ybc_tau[0, 0, :] = t11[ybc_index]
    ybc_tau[0, 1, :] = t12[ybc_index]
    ybc_tau[1, 0, :] = t21[ybc_index]
    ybc_tau[1, 1, :] = t22[ybc_index]

    tau[0, 0, :] = construct_reduced_solution(t11, all_bc_index)
    tau[0, 1, :] = construct_reduced_solution(t12, all_bc_index)
    tau[1, 0, :] = construct_reduced_solution(t21, all_bc_index)
    tau[1, 1, :] = construct_reduced_solution(t22, all_bc_index)
    return


def test_viscous_operator(DxBinv, DyBinv, u, v, sol, nx, ny):

    m = len(u)
    Lap = viscous_operator(DxBinv, DyBinv)
    u = np.append(u, v)
    operator_result = Lap @ u

    # error
    uLinf = la.norm(operator_result[:m] - sol[:m], np.inf)
    vLinf = la.norm(operator_result[m:] - sol[m:], np.inf)

    print(nx, ny, uLinf, vLinf)
    return uLinf, vLinf, 0.0, 0.0


def collocation_ns_momentum(degree, xpts, convection_type, ypts=None, test=None, xboundary='dirichlet'):

    # dimensions
    dim = 2
    # number of derivatives
    deriv = 2
    # runge kutta stage
    stage = 3
    # time
    t = 0.0
    # component
    component = ['u', 'v']

    if test == 'mms':
        # knot vector
        mink = 0.0
        maxk = 2.0 * np.pi
        # length of domain
        L = maxk - mink
        # final time
        T = 0.1
        # time-step
        dt = 5.0e-6

        mu = 1.0
        gamma = 1.4
        M = 1.0
        Re = 1.0
        Pr = 0.72
        function, source, primitive, operator, tau = mms(mu=mu, Re=Re, gamma=gamma, M=M, Pr=Pr, L=L)
    elif test == 'channel':
        # knot vector
        mink = 0.0
        maxk = 2.0
        # length of domain
        L = maxk - mink
        # final time
        T = 30.0
        # time-step
        dt = 1.0e-3

        gamma = 1.4
        M = 1.5
        Re = 3000.0
        Pr = 0.72
        S1 = 110.4
        Tw = 293.15
        Tm = 1.2
        mu = (1.0 + S1/Tw)/(Tm + S1/Tw) * (Tm)**(1.5)

        pressure, temperature = channel_flow(mu=mu, Re=Re, gamma=gamma, M=M, Pr=Pr, L=0.5*L)

    if ypts==None:
        ypts = xpts

    if xboundary == 'periodic':
        xfullu, xders, _ = periodic_init(mink, maxk, xpts, degree, deriv)
        xfullu = xfullu[0:-1]
        xfullB = ss.csc_matrix(periodic_fix(xders[:, 0, :], degree))
        xfullD = ss.csc_matrix(periodic_fix(xders[:, 1, :], degree))
        xfullL = ss.csc_matrix(periodic_fix(xders[:, 2, :], degree))

    elif xboundary == 'dirichlet':
        xfullu, xders, _ = nonperiodic_init(mink, maxk, xpts, degree, deriv)
        xfullB = xders[:, 0, :]
        xfullD = xders[:, 1, :]
        xfullL = xders[:, 2, :]

    yfullu, yders, _ = nonperiodic_init(mink, maxk, ypts, degree, deriv)
    fullu = np.array([xfullu, yfullu])

    yfullB = ss.csc_matrix(yders[:, 0, :])
    yfullD = ss.csc_matrix(yders[:, 1, :])
    yfullL = ss.csc_matrix(yders[:, 2, :])

    # create full 2D mesh
    x, y = np.meshgrid(fullu[0], fullu[1])

    # create full operators by kronecker product
    fullB   = ss.kron(xfullB, yfullB, format='csc')
    fullDx  = ss.kron(xfullD, yfullB, format='csc')
    fullDy  = ss.kron(xfullB, yfullD, format='csc')
    fullDxx = ss.kron(xfullL, yfullB, format='csc')
    fullDyy = ss.kron(xfullB, yfullL, format='csc')
    fullBinv = np.kron(la.inv(xfullB.toarray()), la.inv(yfullB.toarray()))

    yfullLU = ssl.splu(yfullB)
    yfullB = slice_matrix_1d(yfullB)
    xLU = ssl.splu(xfullB)
    yLU = ssl.splu(yfullB)

    del xfullB, yfullB, xfullD, yfullD, xders, yders, xfullu, yfullu, xfullL, yfullL

    DxBinv = fullDx @ fullBinv
    DyBinv = fullDy @ fullBinv
    DxxBinv = fullDxx @ fullBinv
    DyyBinv = fullDyy @ fullBinv

    # able to handle dirichlet/periodic bc_index in x-direction
    _, xbc_index, ybc_index = indexing(fullu, 'bc', xboundary)
    _, xrow_index, xcol_index = indexing(fullu, 'x-nonbc', xboundary)
    _, yrow_index, ycol_index = indexing(fullu, 'y-nonbc', xboundary)
    _, row_index, col_index = indexing(fullu, 'rhs-operator', xboundary)

    # flattened index of all bc points
    if xboundary == 'dirichlet':
        all_bc_index = np.append(xbc_index, ybc_index)

        # bc operators
        xbc_DBinv = np.zeros((dim, len(xrow_index), len(xbc_index)))
        ybc_DBinv = np.zeros((dim, len(yrow_index), len(ybc_index)))

        xbc_DBinv[0, :, :] = slice_matrix(DxBinv, xrow_index, xcol_index)
        ybc_DBinv[0, :, :] = slice_matrix(DxBinv, yrow_index, ycol_index)

        xbc_DBinv[1, :, :] = slice_matrix(DyBinv, xrow_index, xcol_index)
        ybc_DBinv[1, :, :] = slice_matrix(DyBinv, yrow_index, ycol_index)

        xbc_D2Binv = np.zeros((dim, len(xrow_index), len(xbc_index)))
        ybc_D2Binv = np.zeros((dim, len(yrow_index), len(ybc_index)))

        xbc_D2Binv[0, :, :] = slice_matrix(DxxBinv, xrow_index, xcol_index)
        ybc_D2Binv[0, :, :] = slice_matrix(DxxBinv, yrow_index, ycol_index)

        xbc_D2Binv[1, :, :] = slice_matrix(DyyBinv, xrow_index, xcol_index)
        ybc_D2Binv[1, :, :] = slice_matrix(DyyBinv, yrow_index, ycol_index)

        m, n = DxBinv.shape
        assert(m == n)
        mxbc, nxbc = xbc_DBinv[0, :, :].shape
        mybc, nybc = ybc_DBinv[0, :, :].shape
        xbc_LapBinv = np.zeros((dim, dim, mxbc, nxbc))
        ybc_LapBinv = np.zeros((dim, dim, mybc, nybc))

        # i=1, j=1
        xbc_LapBinv[0, 0, :, :] = slice_matrix(LapBinv[:m, :m], xrow_index, xcol_index)
        ybc_LapBinv[0, 0, :, :] = slice_matrix(LapBinv[:m, :m], yrow_index, ycol_index)

        # i=1, j=2
        xbc_LapBinv[0, 1, :, :] = slice_matrix(LapBinv[:m, m:], xrow_index, xcol_index)
        ybc_LapBinv[0, 1, :, :] = slice_matrix(LapBinv[:m, m:], yrow_index, ycol_index)

        # i=2, j=1
        xbc_LapBinv[1, 0, :, :] = slice_matrix(LapBinv[m:, :m], xrow_index, xcol_index)
        ybc_LapBinv[1, 0, :, :] = slice_matrix(LapBinv[m:, :m], yrow_index, ycol_index)

        # i=2, j=2
        xbc_LapBinv[1, 1, :, :] = slice_matrix(LapBinv[m:, m:], xrow_index, xcol_index)
        ybc_LapBinv[1, 1, :, :] = slice_matrix(LapBinv[m:, m:], yrow_index, ycol_index)

    elif xboundary == 'periodic':
        all_bc_index = ybc_index

        # bc operators
        ybc_DBinv = [ss.csc_matrix(slice_matrix(DxBinv, yrow_index, ycol_index)), \
                     ss.csc_matrix(slice_matrix(DyBinv, yrow_index, ycol_index))]

        ybc_D2Binv = [ss.csc_matrix(slice_matrix(DxxBinv, yrow_index, ycol_index)), \
                      ss.csc_matrix(slice_matrix(DyyBinv, yrow_index, ycol_index))]

        ybc_tau = np.zeros((dim, dim, len(ybc_index)))

    # operators
    D = [ss.csc_matrix(slice_matrix(fullDx, row_index, col_index)), \
         ss.csc_matrix(slice_matrix(fullDy, row_index, col_index))]
    D2 = [ss.csc_matrix(slice_matrix(fullDxx, row_index, col_index)), \
          ss.csc_matrix(slice_matrix(fullDyy, row_index, col_index))]

    del fullB, fullDxx, fullDyy, fullBinv, DxBinv, DyBinv, DxxBinv, DyyBinv #, fullDx, fullDy

    # time-step vectors
    rka, rkb, rkc = low_storage_rk(stage)

    nx = len(fullu[0])
    ny = len(fullu[1])
    N = nx*(ny-2)
    NN = nx*ny

    # registers for low-storage rk
    U2 = np.zeros((dim, N))
    U1 = np.zeros_like(U2)
    tempU = np.zeros_like(U2)
    tempfullU = np.zeros((dim, NN))
    Rho2 = np.zeros(N)
    Rho1 = np.zeros_like(Rho2)
    tempRho = np.zeros_like(Rho2)
    T2 = np.zeros(N)
    T1 = np.zeros_like(T2)
    tempT = np.zeros_like(T2)
    S = np.zeros((4, N))
    sol = np.zeros((dim, N))
    tau = np.zeros((dim, dim, N))
    Linf = np.zeros(dim)
    force = np.zeros((dim, N))

    # initialize alpha
    alpha = np.zeros((dim, N))
    alpha_T = np.zeros(N)
    alpha_rho_u = np.zeros(N)
    alpha_rho_v = np.zeros(N)
    alpha_rho_T = np.zeros(N)
    alpha_conv_T = np.zeros((dim, N))
    alpha_conv_div = np.zeros((dim, dim, N))
    alpha_tau = np.zeros((dim, dim, N))
    alpha_full_u = np.zeros((dim, NN))

    # boundary condition from operators
    cont_bc = np.zeros(N)
    conv_bc = np.zeros((dim, N))
    Lap_bc = np.zeros((dim, N))
    rho_T_bc = np.zeros((dim, N))
    Tconv_bc = np.zeros(N)
    TLap_bc =  np.zeros(N)
    Tconv_bc = np.zeros(N)
    Tdivu_bc = np.zeros(N)
    Tdudx_bc = np.zeros((dim, dim, N))

    # boundary condition vectors
    u_n_ybc = np.zeros((dim, len(ybc_index)))
    T_ybc = np.zeros(len(ybc_index))
    rho_T_ybc = np.zeros(len(ybc_index))
    rho_u_ybc = np.zeros(len(ybc_index))
    rho_v_ybc = np.zeros(len(ybc_index))
    conv_div_ybc = np.zeros((dim, dim, len(ybc_index)))
    rho_T_u_ybc = np.zeros(len(ybc_index))
    rho_T_v_ybc = np.zeros(len(ybc_index))

    ### test operators
    if test == 'viscous' and xboundary == 'periodic':
        for d in range(2):
            tempU[d, :] = construct_reduced_solution(function[d+1](x, y, t).flatten('F'), all_bc_index)
            u_n_ybc[d, :] = yboundary(function[d+1], t, fullu, ybc_index)

        x_periodic_velocity_gradient_bc(tau, ybc_tau,
                tempU, u_n_ybc, fullu, xboundary, xLU, yfullLU, tempfullU,
                all_bc_index, ybc_index, fullDx, fullDy, alpha_full_u, mu)

        for d in range(2):
            for e in range(2):
                alpha_tau[d, e, :] = LU.solve(tau[d, e, :])

        # dti1/dx + dti2/dy
        viscous_bc_u = ybc_DBinv[0, :, :] @ ybc_tau[0, 0, :] + ybc_DBinv[1, :, :] @ ybc_tau[0, 1, :]
        viscous_bc_v = ybc_DBinv[0, :, :] @ ybc_tau[1, 0, :] + ybc_DBinv[1, :, :] @ ybc_tau[1, 1, :]

        viscous_u = D[0] @ alpha_tau[0, 0, :] + D[1] @ alpha_tau[0, 1, :] + viscous_bc_u
        viscous_v = D[0] @ alpha_tau[1, 0, :] + D[1] @ alpha_tau[1, 1, :] + viscous_bc_v

        solu = construct_reduced_solution(operator[0](x, y, t).flatten('F'), all_bc_index)
        solv = construct_reduced_solution(operator[1](x, y, t).flatten('F'), all_bc_index)

        error_u = la.norm(viscous_u - solu, np.inf)
        error_v = la.norm(viscous_v - solv, np.inf)

        print(xpts, error_u, error_v)
        return error_u, error_v, 0.0, 0.0

    elif test == 'convection':
        convection = convection_eval(DxBinv, DyBinv, rho, u, v, convection_type)
        return test_convection_eval(convection, convection_type, x, y, t)

    if test == 'mms':
        rho_ybc = None

        for d in range(dim):
            # initialize u
            U2[d, :] = construct_reduced_solution(function[d+1](x, y, t).flatten('F'), all_bc_index)
        Rho2[:] = construct_reduced_solution(function[0](x, y, t).flatten('F'), all_bc_index)
        T2[:] = construct_reduced_solution(function[3](x, y, t).flatten('F'), all_bc_index)

    elif test == 'channel':
        function = None
        primitive = None

        T2[:] = Tm
        U2[0, :] = 1.0 # M * np.sqrt(gamma * 286.0 * 293.15)
        U2[1, :] = 0.0
        Pstar = pressure(x, y, t).flatten('F') / np.mean(pressure(x, y, t).flatten('F'))
        fullrho = Pstar * gamma * M**2 / np.mean(T2)
        Rho2[:] = construct_reduced_solution(fullrho, all_bc_index)
        rho_ybc = fullrho[ybc_index]
        sum_rho0 = np.sum(fullrho)

    # number of time-steps
    N = int(np.rint(T / dt))

    for time in range(N):

        # SSP Runge-Kutta
        for m in range(1, stage+1):
            rkt = t + rkc[m]*dt

            Rho1[:] = Rho2 + dt * (rka[m, m-1] - rkb[m-1]) * Rho1
            U1[0, :] = U2[0, :] + dt * (rka[m, m-1] - rkb[m-1]) * U1[0, :]
            U1[1, :] = U2[1, :] + dt * (rka[m, m-1] - rkb[m-1]) * U1[1, :]
            T1[:] = T2 + dt * (rka[m, m-1] - rkb[m-1]) * T1

            tempU[:, :] = U1.copy()
            tempT[:] = T1.copy()
            tempRho[:] = Rho1.copy()

            ### FORM BC
            periodic_boundary_condition(u_n_ybc, T_ybc, rho_T_ybc, rho_u_ybc, rho_v_ybc,
                                conv_div_ybc, rho_T_u_ybc, rho_T_v_ybc,
                                rkt, fullu, ybc_index, function, primitive, test, rho_ybc=rho_ybc)
            # continuity bc
            continuity_bc(cont_bc, ybc_DBinv, rho_u_ybc, rho_v_ybc)

            # momentum bc
            momentum_bc(conv_bc, Lap_bc, rho_T_bc, tau, ybc_tau,
                fullu, xLU, yfullLU, tempfullU, xboundary,
                tempU, tempRho, u_n_ybc, rho_T_ybc, conv_div_ybc, rho_u_ybc, rho_v_ybc,
                all_bc_index, ybc_index, fullDx, fullDy, alpha_full_u, mu,
                ybc_DBinv)

            # energy bc
            energy_bc(Tconv_bc, TLap_bc, Tdudx_bc, Tdivu_bc,
                    tempU, tempT, tempRho,
                    u_n_ybc, T_ybc, rho_u_ybc, rho_v_ybc, rho_T_u_ybc, rho_T_v_ybc,
                    ybc_DBinv, ybc_D2Binv)

            ### Solve linear system
            linear_system(xLU, yLU, Rho1, U1[0, :], U1[1, :], T1, tau,
                          alpha, alpha_T, alpha_tau,
                          alpha_rho_u, alpha_rho_v, alpha_rho_T,
                          alpha_conv_T, alpha_conv_div)

            ### CONTINUITY EQUATION
            # initialize source term
            if test == 'mms':
                S[0, :] = construct_reduced_solution(source[0](x, y, rkt).flatten('F'), all_bc_index)

            # (m-stage, m order) Low-storage Explict Runge-Kutta
            Rho1[:] = - (D[0] @ alpha_rho_u +  \
                         D[1] @ alpha_rho_v + cont_bc) + \
                         S[0, :]

            ### MOMENTUM EQUATION
            for d in range(dim):
                # initialize source term
                if test == 'mms':
                    S[d+1, :] = construct_reduced_solution(source[d+1](x, y, rkt).flatten('F'), all_bc_index)

                elif test == 'channel':
                    force[0, :] = - (np.mean(np.abs(ybc_tau[0, 1, :])) / Re)
                    force[1, :] = 0.0

                # (m-stage, m order) Low-storage Explict Runge-Kutta
                U1[d, :] = (-0.5/tempRho) *    \
                           (tempRho * (tempU[0, :] * (D[0] @ alpha[d, :]) +  \
                            tempU[1, :] * (D[1] @ alpha[d, :])) + \
                            D[0] @ alpha_conv_div[d, 0, :] + D[1] @ alpha_conv_div[d, 1, :] - \
                            tempU[d, :] * (D[0] @ alpha_rho_u + D[1] @ alpha_rho_v) + conv_bc[d, :]) +  \
                           (1.0/(Re * tempRho)) *                        \
                           (D[0] @ alpha_tau[d, 0, :] +  \
                            D[1] @ alpha_tau[d, 1, :] + Lap_bc[d, :]) - \
                           (1.0 / (gamma * M**2 * tempRho)) * \
                           (D[d] @ alpha_rho_T + rho_T_bc[d, :]) + \
                            S[d+1, :] + force[d, :]

            ### ENERGY EQUATION
            # initialize source term
            if test == 'mms':
                S[3, :] = construct_reduced_solution(source[3](x, y, rkt).flatten('F'), all_bc_index)

            # (m-stage, m order) Low-storage Explict Runge-Kutta
            T1[:] = (-0.5/tempRho) *    \
                    (tempRho * (tempU[0, :] * (D[0] @ alpha_T) +  \
                                tempU[1, :] * (D[1] @ alpha_T)) + \
                     D[0] @ alpha_conv_T[0, :] + D[1] @ alpha_conv_T[1, :] - \
                     tempT * (D[0] @ alpha_rho_u + D[1] @ alpha_rho_v) + Tconv_bc) - \
                    (gamma - 1.0) * (tempT * (D[0] @ alpha[0, :] + D[1] @ alpha[1, :]) + Tdivu_bc) + \
                    (gamma * (gamma - 1.0) * M**2 / (Re * tempRho)) * \
                    (tau[0, 0, :] * (D[0] @ alpha[0, :] + Tdudx_bc[0, 0, :]) + \
                     tau[0, 1, :] * (D[1] @ alpha[0, :] + Tdudx_bc[0, 1, :]) + \
                     tau[1, 0, :] * (D[0] @ alpha[1, :] + Tdudx_bc[1, 0, :]) + \
                     tau[1, 1, :] * (D[1] @ alpha[1, :] + Tdudx_bc[1, 1, :])) + \
                    (gamma * mu) / (Re * Pr * tempRho) * \
                    (D2[0] @ alpha_T + D2[1] @ alpha_T + TLap_bc) + \
                     S[3, :]

            Rho2[:] = Rho2 + rkb[m] * dt * Rho1
            U2[0, :] = U2[0, :] + rkb[m] * dt * U1[0, :]
            U2[1, :] = U2[1, :] + rkb[m] * dt * U1[1, :]
            T2[:] = T2 + rkb[m] * dt * T1

        t = (time + 1) * dt

        if test == 'mms':
            for d in range(dim):
                # exact solution
                sol[d, :] = construct_reduced_solution(function[d+1](x, y, t).flatten('F'), all_bc_index)

                # L_inf error
                Linf[d] = la.norm(sol[d, :] - U2[d, :], np.inf)
                assert(Linf[d] < 1e2)

            solT = construct_reduced_solution(function[3](x, y, t).flatten('F'), all_bc_index)
            LinfT = la.norm(solT - T2, np.inf)

            solrho = construct_reduced_solution(function[0](x, y, t).flatten('F'), all_bc_index)
            Linfrho = la.norm(solrho - Rho2, np.inf)

            print('Time = %f, rhoLinf = %e, uLinf = %e, vLinf = %e, TLinf = %e, Points = %d x %d' % (t, Linfrho, Linf[0], Linf[1], LinfT, xpts, ypts))

        elif test == 'channel':
            fullrho = construct_full_solution(Rho2, None, rho_ybc, fullu, xboundary)
            fullU = construct_full_solution(U2[0, :], None, u_n_ybc[0, :], fullu, xboundary)
            fullV = construct_full_solution(U2[1, :], None, u_n_ybc[1, :], fullu, xboundary)
            fullT = construct_full_solution(T2, None, T_ybc, fullu, xboundary)

            print('Time = %f, <rho> = %e, <u> = %e, <v> = %e, <T> = %e, Points = %d x %d'
                    % (t, np.mean(fullrho), np.mean(fullU), np.mean(fullV), np.mean(fullT), xpts, ypts))

    if test == 'mms':
        return Linfrho, Linf[0], Linf[1], LinfT
    elif test == 'channel':
        xfullu, _, _ = periodic_init(mink, maxk, xpts, degree, 1)
        yfullu, _, _ = nonperiodic_init(mink, maxk, ypts, degree, 1)
        return xfullu, yfullu, fullrho, fullU, fullV, fullT


def MMS_test_collocation_ns_2D(degree):
    import csv

    test = 'mms'

    error1 = []
    error2 = []
    error3 = []
    error4 = []
    xpts = np.array([40, 50, 60, 70, 80, 90, 100, 110, 120])

    for i in range(len(xpts)):

        a, b, c, d = collocation_ns_momentum(degree, xpts[i], 'convection skew-symmetric', xboundary='periodic', test=test)

        error1.append(a)
        error2.append(b)
        error3.append(c)
        error4.append(d)

        temp = zip(xpts, error1, error2, error3, error4)

        with open('collocation-ns-momentum-rk3-degree%d.csv'% degree, "w") as output:
            writer = csv.writer(output, lineterminator='\n', delimiter=',')
            for val in temp:
                writer.writerow(val)


def CHANNEL_test_collocation_ns_2D(degree):
    import csv

    test = 'channel'

    xpts = np.array([30])

    for i in range(len(xpts)):

        x, y, a, b, c, d = collocation_ns_momentum(degree, xpts[i], 'convection skew-symmetric', xboundary='periodic', test=test)

        temp = zip(a, b, c, d)

        with open('RE150-collocation-ns-momentum-rk3-degree%d.csv'% degree, "w") as output:
            writer = csv.writer(output, lineterminator='\n', delimiter=',')
            for val in temp:
                writer.writerow(val)

        temp = zip(x, y)
        with open('mesh-degree%d.csv'% degree, "w") as output:
            writer = csv.writer(output, lineterminator='\n', delimiter=',')
            for val in temp:
                writer.writerow(val)


def post_process():
    import csv
    import string
    import itertools
    import scipy.stats.mstats as stats


    x = []; y = []; rho = []; u = []; v = []; T = []
    f = []; m = []

    label = ['\\rho', 'u', 'v', 'T']
    Re = 150
    degree = 7

    f.append('RE%d-collocation-ns-momentum-rk3-degree%d.csv'% (Re, degree))
    m.append('mesh-degree%d.csv'% degree)

    with open(f[0], 'r') as data:
        reader = csv.reader(data, delimiter=',')
        for row in reader:
            rho.append(np.float(row[0]))
            u.append(np.float(row[1]))
            v.append(np.float(row[2]))
            T.append(np.float(row[3]))

    with open(m[0], 'r') as data:
        reader = csv.reader(data, delimiter=',')
        for row in reader:
            x.append(np.float(row[0]))
            y.append(np.float(row[1]))

    #rho = np.reshape(rho, (len(x), len(y)), 'F')
    #u = np.reshape(u, (len(x), len(y)), 'F')
    #v = np.reshape(v, (len(x), len(y)), 'F')
    #T = np.reshape(T, (len(x), len(y)), 'F')

    xboundary = 'periodic'
    xpts = len(x)
    ypts = xpts
    mink = np.min(x)
    maxk = np.max(x)
    deriv = 1

    if xboundary == 'periodic':
        xfullu, xders, _ = periodic_init(mink, maxk, xpts, degree, deriv)
        xfullB = periodic_fix(xders[:, 0, :], degree)
        xfullD = periodic_fix(xders[:, 1, :], degree)

    elif xboundary == 'dirichlet':
        xfullu, xders, _ = nonperiodic_init(mink, maxk, xpts, degree, deriv)
        xfullB = xders[:, 0, :]
        xfullD = xders[:, 1, :]

    yfullu, yders, _ = nonperiodic_init(mink, maxk, ypts, degree, deriv)
    fullu = np.array([xfullu, yfullu])

    yfullB = yders[:, 0, :]
    yfullD = yders[:, 1, :]

    # create full 2D mesh
    xx, yy = np.meshgrid(fullu[0], fullu[1])

    # create full operators by kronecker product
    fullB = np.kron(xfullB, yfullB)
    fullDx = np.kron(xfullD, yfullB)
    fullDy = np.kron(xfullB, yfullD)

    vorticity = fullDx @ v[:-30] - fullDy @ u[:-30]
    vorticity = np.reshape(np.append(vorticity, vorticity[:30]), (len(xfullu), len(yfullu)), 'F')

    plt.quiver(xx, yy, vorticity, 0.0)
    #plt.pcolormesh(xx, yy, vorticity)
    plt.show()


def main():
    # collocation_ns_momentum(5, 20)
    #CHANNEL_test_collocation_ns_2D(3)
    MMS_test_collocation_ns_2D(7)
    #post_process()


if __name__ == "__main__":
    main()
