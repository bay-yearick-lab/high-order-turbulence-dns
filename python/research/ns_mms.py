import numpy as np
import sympy as sp


def channel_flow(w=1.0, u0=0.05, v0=0.08, r0=0.1, T0=0.5, eps=0.5, mu=1.0, Re=1.0, gamma=1.4, M=1.0, Pr=0.72, L=1.0):

    x, y, t = sp.symbols('x, y, t')

    dpdx = 10.0
    dTdx = -dpdx

    P = - dpdx * x + (dpdx * L + 100.5)
    T = - dTdx * x + (dTdx * L + 293.15)

    P = sp.lambdify((x, y, t), P, 'numpy')
    T = sp.lambdify((x, y, t), T, 'numpy')

    return P, T


def mms(w=1.0, u0=0.05, v0=0.08, r0=0.1, T0=0.5, eps=0.5, mu=1.0, Re=1.0, gamma=1.4, M=1.0, Pr=0.72, L=1.0):

    x, y, t = sp.symbols('x, y, t')

    u = u0 * (sp.cos(4.0 * np.pi * x / L) * sp.cos(2.0 * np.pi * y / L) * sp.exp(-t**2))
    v = v0 * (sp.cos(4.0 * np.pi * x / L) * sp.cos(2.0 * np.pi * y / L) * sp.exp(-t**2))
    rho = r0 * (sp.cos(4.0 * np.pi * x / L) * sp.sin(2.0 * np.pi * y / L) * sp.exp(-t**2) + 1.5)
    T = T0 * (sp.sin(4.0 * np.pi * x / L) * sp.cos(2.0 * np.pi * y / L) * sp.exp(-t**2) + 1.5)

    rhou = rho * u
    rhov = rho * v
    rhoT = rho * T
    rhouu = rho * u**2
    rhouv = rho * u * v
    rhovv = rho * v**2
    rhoTu = rho * T * u
    rhoTv = rho * T * v

    div_u = (rho * u**2).diff(x) + (rho * u * v).diff(y)

    div_v = (rho * u * v).diff(x) + (rho * v**2).diff(y)

    adv_u = rho * u * u.diff(x) + rho * v * u.diff(y)

    adv_v = rho * u * v.diff(x) + rho * v * v.diff(y)

    skew_u = 0.5 * (rho * u * u.diff(x) + rho * v * u.diff(y) + \
                   (rho * u**2).diff(x) + (rho * u * v).diff(y) - \
                   u * ((rho * u).diff(x) + (rho * v).diff(y)))

    skew_v = 0.5 * (rho * u * v.diff(x) + rho * v * v.diff(y) + \
                   (rho * u * v).diff(x) + (rho * v**2).diff(y) - \
                   v * ((rho * u).diff(x) + (rho * v).diff(y)))

    skew_T = 0.5 * (rho * u * T.diff(x) + rho * v * T.diff(y) + \
                   (rho * T * u).diff(x) + (rho * T * v).diff(y) - \
                   T * ((rho * u).diff(x) + (rho * v).diff(y)))

    t11 = mu * (4.0/3.0 * u.diff(x) - 2.0/3.0 * v.diff(y))
    t12 = mu * (u.diff(y) + v.diff(x))
    t21 = t12
    t22 = mu * (4.0/3.0 * v.diff(y) - 2.0/3.0 * u.diff(x))

    viscous_u = (t11).diff(x) + (t12).diff(y)
    viscous_v = (t21).diff(x) + (t22).diff(y)

    S1 = rho.diff(t) + (rho * u).diff(x) + (rho * v).diff(y)
    S2u = u.diff(t) + (1.0/rho) * skew_u - 1.0/(Re * rho) * viscous_u + 1.0/(gamma * M**2 * rho) * (rho * T).diff(x)
    S2v = v.diff(t) + (1.0/rho) * skew_v - 1.0/(Re * rho) * viscous_v + 1.0/(gamma * M**2 * rho) * (rho * T).diff(y)
    S3 = T.diff(t) + (1.0/rho) * skew_T + (gamma - 1.0) * T * (u.diff(x) + v.diff(y)) - \
         gamma * (gamma - 1.0) * M**2 / (Re * rho) * (t11 * u.diff(x) + t12 * u.diff(y) + t21 * v.diff(x) + t22 * v.diff(y)) - \
         gamma * mu / (Re * Pr * rho) * (T.diff(x, x) + T.diff(y, y))

    rho = sp.lambdify((x, y, t), rho, 'numpy')
    u = sp.lambdify((x, y, t), u, 'numpy')
    v = sp.lambdify((x, y, t), v, 'numpy')
    T = sp.lambdify((x, y, t), T, 'numpy')

    S1 = sp.lambdify((x, y, t), S1, 'numpy')
    S2u = sp.lambdify((x, y, t), S2u, 'numpy')
    S2v = sp.lambdify((x, y, t), S2v, 'numpy')
    S3 = sp.lambdify((x, y, t), S3, 'numpy')

    rhou = sp.lambdify((x, y, t), rhou, 'numpy')
    rhov = sp.lambdify((x, y, t), rhov, 'numpy')
    rhouu = sp.lambdify((x, y, t), rhouu, 'numpy')
    rhouv = sp.lambdify((x, y, t), rhouv, 'numpy')
    rhovv = sp.lambdify((x, y, t), rhovv, 'numpy')
    rhoT = sp.lambdify((x, y, t), rhoT, 'numpy')
    rhoTu = sp.lambdify((x, y, t), rhoTu, 'numpy')
    rhoTv = sp.lambdify((x, y, t), rhoTv, 'numpy')

    viscous_u = sp.lambdify((x, y, t), viscous_u, 'numpy')
    viscous_v = sp.lambdify((x, y, t), viscous_v, 'numpy')

    t11 = sp.lambdify((x, y, t), t11, 'numpy')
    t12 = sp.lambdify((x, y, t), t12, 'numpy')
    t21 = sp.lambdify((x, y, t), t21, 'numpy')
    t22 = sp.lambdify((x, y, t), t22, 'numpy')

    return [rho, u, v, T], [S1, S2u, S2v, S3], [rhou, rhov, rhouu, rhouv, rhovv, rhoT, rhoTu, rhoTv], [viscous_u, viscous_v], [t11, t12, t21, t22]
