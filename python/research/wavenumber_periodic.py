import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import matplotlib.pyplot as plt
from findspan import findspan
from derivatives import derivatives
from collocation import collocation_points
from knot_multiplicity import periodic_knot_multiplicity


def BB(u, a, x, n0, n1, n_splines, U, degree, band):

    n = max(n0, n1)
    i = findspan(n_splines - 1, degree, u, U)
    derivatives(i, u, degree, n, U, x)

    return x[n0, a] * x[n1, a + band]


def periodic_fix(M, degree, oversample=1):

    m, n = M.shape
    if (m == (n - degree + 1)):
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
            row = np.zeros((degree - 1, m))

            for i in range(degree):
                M = np.delete(M, -1, 1)
                if (i < (degree-1)):
                    n = degree - (i + 1)
                    row[i, :] = M[-(i+1), n:]
                    row[i, :] = np.roll(row[i, :], n)

            for i in range(1, degree):
                M[-i, :] = row[i-1, :]
    else:
        if (degree%2 == 0):
            print("NOT IMPLEMENTED YET")
            quit()

        else:
            m, n = M.shape
            row = np.zeros((degree - 1, m))

            for i in range(degree - 1):
                M = np.delete(M, -1, 1)
        M[-degree:, :] = np.flip(M[:degree, :], 0)

    return M


def compact_fd_penta(a, N):

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    U = np.linspace(mink, maxk, N)
    # points
    pts = N - 1
    dx = U[1] - U[0]

    bw = 3
    m = np.zeros(2 * bw + 1)
    d = np.zeros(2 * bw + 1)

    d[bw-3] = -1.0/50.0 * (9.0 * a - 4.0) / (6.0 * dx)
    d[bw-2] = -1.0/150.0 * (568.0 * a - 183.0) / (4.0 * dx)
    d[bw-1] = -1.0/6.0 * (12.0 - 7.0 * a) / (2.0 * dx)
    d[bw] = 0.0
    d[bw+1] = -d[bw-1]
    d[bw+2] = -d[bw-2]
    d[bw+3] = -d[bw-3]

    m[bw-2] = 1.0/20.0 * (-3.0 + 8.0 * a)
    m[bw-1] = a
    m[bw] = 1.0
    m[bw+1] = a
    m[bw+2] = m[bw-2]

    kmax = np.pi / dx
    kexact = np.linspace(0.0, np.pi, pts)

    numodd = 0.0
    den = m[bw]
    for i in range(1, bw + 1):
        numodd += 2.0 * d[bw+i] * np.sin(kexact * i)
        den += 2.0 * m[bw+i] * np.cos(kexact * i)

    lam = abs(numodd / den)

#    D = np.zeros((pts, pts))
#    row = np.roll(d, -4)
#    print(row)
#    for i in range(pts):
#        row = np.roll(row, 1)
#        D[i, :] = row
#    print(D)

#    plt.figure(1, figsize=(14, 12))
#    plt.plot(kexact/max(kexact), lam/kmax, '--', linewidth=3, label='6th-Order Compact')
#    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
#    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=36)
#    plt.xlabel('$k/k_{max}$', fontsize=36)
#    plt.xticks(fontsize=36)
#    plt.yticks(fontsize=36)
#    plt.legend(loc='best', fontsize=36)
#    plt.tight_layout()
#    plt.show()

    return kexact, lam, kmax


def compact_fd_tri(a, N):

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    U = np.linspace(mink, maxk, N)
    # points
    pts = N - 1
    dx = U[1] - U[0]

    bw = 3
    m = np.zeros(2 * bw + 1)
    d = np.zeros(2 * bw + 1)

    d[bw-2] = -1.0/3.0 * (4.0 * a - 1.0) / (4.0 * dx)
    d[bw-1] = -2.0/3.0 * (a + 2.0) / (2.0 * dx)
    d[bw] = 0.0
    d[bw+1] = -d[bw-1]
    d[bw+2] = -d[bw-2]

    m[bw-1] = a
    m[bw] = 1.0
    m[bw+1] = a


#    d[bw-3] = -1.0/10.0 * (-3.0 * a + 1.0) / (6.0 * dx)
#    d[bw-2] = -1.0/15.0 * (32.0 * a - 9.0) / (4.0 * dx)
#    d[bw-1] = -1.0/6.0 * (a + 9.0) / (2.0 * dx)
#    d[bw] = 0.0
#    d[bw+1] = -d[bw-1]
#    d[bw+2] = -d[bw-2]
#    d[bw+3] = -d[bw-3]
#
#    m[bw-1] = a
#    m[bw] = 1.0
#    m[bw+1] = a

    kmax = np.pi / dx
    kexact = np.linspace(0.0, np.pi, pts)

    numodd = 0.0
    den = m[bw]
    for i in range(1, bw + 1):
        numodd += 2.0 * d[bw+i] * np.sin(kexact * i)
        den += 2.0 * m[bw+i] * np.cos(kexact * i)

    lam = abs(numodd / den)

#    D = np.zeros((pts, pts))
#    row = np.roll(d, -4)
#    print(row)
#    for i in range(pts):
#        row = np.roll(row, 1)
#        D[i, :] = row
#    print(D)

#    plt.figure(1, figsize=(14, 12))
#    plt.plot(kexact/max(kexact), lam/kmax, '--', linewidth=3, label='6th-Order Compact')
#    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
#    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=36)
#    plt.xlabel('$k/k_{max}$', fontsize=36)
#    plt.xticks(fontsize=36)
#    plt.yticks(fontsize=36)
#    plt.legend(loc='best', fontsize=36)
#    plt.tight_layout()
#    plt.show()

    return kexact, lam, kmax


def galerkin_periodic(deriv, degree, N):
    from scipy.integrate import quad
    from knot_multiplicity import knot_multiplicity

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    U = np.linspace(mink, maxk, N)
    dx = U[1] - U[0]
    # knot multiplicity
    u = knot_multiplicity(degree, U)
    U = periodic_knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # points
    pts = N - 1

    x = np.zeros((deriv + 1, degree + 1))

    # parameter
    zero = 0

    M = np.zeros((pts, pts))
    D = np.zeros((pts, pts))

    # position of B_i in respective subintervals
    a = np.arange(degree + 1)[::-1]

    xi = np.zeros(degree + 2)

    # quad rule with B-splines
    for i in range(pts):

        # coordinates of subintervals
        for j in range(len(xi)):
            xi[j] = u[i + j]

        for band in range(-degree, degree + 1):
            for j in range(len(a) - abs(band)):

                if (band >= 0 and (i + band) < pts):
                    k = j + band

                    # (B_j, B_i')
                    D[i, i + band] += quad(BB, xi[k], xi[k+1], args=(a[k], x, zero, deriv, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    M[i, i + band] += quad(BB, xi[k], xi[k+1], args=(a[k], x, zero, zero, n_splines, U, degree, band))[0]

                elif (band < 0 and (i + band) >= 0):
                    # (B_j, B_i')
                    D[i, i + band] += quad(BB, xi[j], xi[j+1], args=(a[j], x, zero, deriv, n_splines, U, degree, band))[0]

                    # (B_j, B_i)
                    M[i, i + band] += quad(BB, xi[j], xi[j+1], args=(a[j], x, zero, zero, n_splines, U, degree, band))[0]

    m, n = M.shape
    assert(m == n)

    m = np.zeros((degree + 1, n))
    d = np.zeros((degree + 1, n))

#    kmax = max(kexact)
    kmax = np.pi/dx
    kexact = np.linspace(0.0, np.pi, n)
    mrow= np.roll(M[degree, :(2 * degree + 1)], -degree)
    drow = np.roll(D[degree, :(2 * degree + 1)], -degree)


    for i in range(degree + 1):
        m[i, :] = mrow[i]
        d[i, :] = drow[i]

    numodd = 0.0
    den = m[0, :]
    numeven = d[0, :]
    for i in range(1, degree + 1):
        numodd += 2.0 * d[i, :] * np.sin(kexact * i)
        den += 2.0 * m[i, :] * np.cos(kexact * i)

        numeven += 2.0 * d[i, :] * np.cos(kexact * i)

    if (deriv == 1):
        lam = abs(numodd / den)
    elif (deriv == 2):
        lam = np.sqrt(-numeven / den)

#    lam = 2.0 * np.pi * lam / N

#    plt.figure(1, figsize=(14, 12))
#    plt.plot(kexact/kmax, lam/kmax, linewidth=3, label='$B_i^{%d}$'%degree)
#    plt.plot(kexact/kmax, kexact/kmax,'k', linewidth=4, label='Exact')
#    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=36)
#    plt.xlabel('$k/k_{max}$', fontsize=36)
#    plt.xticks(fontsize=36)
#    plt.yticks(fontsize=36)
#    plt.legend(loc='best', fontsize=36)
#    plt.tight_layout()
#    plt.show()

    return kexact, lam, kmax

def V(lamda, M, P, I, u):
    A = M @ (la.inv(M.T @ M + lamda * P) @ M.T)
    V = la.norm(u -  A @ u, 2)**2
    V = V / np.trace(I - A)**2
    return V*u.shape[0]

def collocation_smooth_periodic(lamda, degree, N, type, p_deriv, random):
    from scipy.integrate import quad
    from scipy.optimize import minimize_scalar
    from knot_multiplicity import knot_multiplicity

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    deriv = 1
    U = np.linspace(mink, maxk, N)
    # knot multiplicity
    U = periodic_knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = p_deriv
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

    fullM = ders[:, 0, :]
    fullD = ders[:, 1, :]

    M = periodic_fix(fullM, degree)
    D = periodic_fix(fullD, degree)
    dx = u[1] - u[0]

    # Galerkin
    U = np.linspace(mink, maxk, N+1)
    dx = U[1] - U[0]
    # knot multiplicity
    u = knot_multiplicity(degree, U)
    U = periodic_knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # points
    pts = N-1

    x = np.zeros((n + 1, degree + 1))

    # parameter
    zero = 0

    P = np.zeros((pts, pts))

    # position of B_i in respective subintervals
    a = np.arange(degree + 1)[::-1]

    xi = np.zeros(degree + 2)

    # quad rule with B-splines
    for i in range(pts):

        # coordinates of subintervals
        for j in range(len(xi)):
            xi[j] = u[i + j]

        for band in range(-degree, degree + 1):
            for j in range(len(a) - abs(band)):

                if (band >= 0 and (i + band) < pts):
                    k = j + band

                    # (B_j'', B_i'')
                    P[i, i + band] += quad(BB, xi[k], xi[k+1], args=(a[k], x, n, n, n_splines, U, degree, band))[0]

                elif (band < 0 and (i + band) >= 0):
                    # (B_j'', B_i'')
                    P[i, i + band] += quad(BB, xi[j], xi[j+1], args=(a[j], x, n, n, n_splines, U, degree, band))[0]

    if (degree%2 == 0):
        bw = int(np.ceil((degree - 1) / 2))
    else:
        bw = int((degree - 1) / 2)

    kexact = np.linspace(0.0, np.pi, pts)
    x = np.linspace(0.0, 2*np.pi, pts)
    lam = np.zeros(pts)

    D = np.roll(D, -bw, 1)
    M = np.roll(M, -bw, 1)

    Prow = P[degree, :]

    for i in range(P.shape[0]):
        P[i, :] = np.roll(Prow, -degree + i)

    if lamda == 'original':
        D = D @ la.inv(M)
        for i in range(pts):
            lam[i] = np.median(np.abs(D @ np.exp(i * 1j * x/2)))
        return kexact, lam, N/2

    D = D @ (la.inv(M.T @ M + lamda/N * P) @ M.T)

    if type == 'wavenumber':
        for i in range(pts):
            lam[i] = np.median(np.abs(D @ np.exp(i * 1j * x/2)))

        #plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
        #plt.plot(kexact/max(kexact), lam/(N/2), linewidth=3, label='$B_i^{%d}$'%degree)
        #plt.show()

        return kexact, lam, N/2
    elif type == 'objective':
        param = np.linspace(0.0, 0.1, 10000)
        u = np.sin(co[0, 1:])
        for i in range(len(co[0,1:-1])):
            u[i] = u[i] + random[i]

        I = np.eye(u.shape[0])
        objective = np.zeros(len(param))

        for i in range(len(param)):
            objective[i] = V(param[i], M, P, I, u)

        res = minimize_scalar(V, args=(M, P, I, u), bounds=(0.0,1.0), method='bounded')
        op = res.x
        func = res.fun

        #plt.plot(co[0,1:], u)
        #plt.plot(co[0,1:], M @ (la.inv(M.T @ M + op * P) @ M.T)@u)
        #plt.show()
        #print(op)

        #plt.ylim(top=np.max(objective[1:]*1.2))
        #plt.xlim(left=0.0, right=0.1)
        #plt.plot(param, objective)
        #plt.plot(op, func, 'ro')
        #plt.show()

        return co[0,1:], M @ (la.inv(M.T @ M + op * P) @ M.T)@u, param, objective, op, func


def collocation_periodic(deriv, degree, N):

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    U = np.linspace(mink, maxk, N)
    # knot multiplicity
    U = periodic_knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = 2
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

    fullM = ders[:, 0, :]
    fullD = ders[:, 1, :]
    fullL = ders[:, 2, :]

    M = periodic_fix(fullM, degree)
    D = periodic_fix(fullD, degree)
    L = periodic_fix(fullL, degree)
    dx = u[1] - u[0]

    if (degree%2 == 0):
        bw = int(np.ceil((degree - 1) / 2))
    else:
        bw = int((degree - 1) / 2)

    M = np.roll(M, -bw, 1)
    D = np.roll(D, -bw, 1)
    L = np.roll(L, -bw, 1)

    m, n = M.shape
    assert(m == n)

    m = np.zeros((bw + 1, n))
    d = np.zeros((bw + 1, n))
    l = np.zeros((bw + 1, n))

    for i in range(bw + 1):
        m[i, :] = M[0, i]
        d[i, :] = D[0, i]
        l[i, :] = L[0, i]

    kexact = np.linspace(0.0, np.pi, n)
#    kmax = max(kexact)

    if (degree%2 == 0):
        kmax = np.pi/(2.0 * dx)
    else:
        kmax = np.pi/dx


    numodd = 0.0
    den = m[0, :]
    numeven = l[0, :]
    for i in range(1, bw + 1):
        numodd += 2.0 * d[i, :] * np.sin(kexact * i)
        den += 2.0 * m[i, :] * np.cos(kexact * i)

        numeven += 2.0 * l[i, :] * np.cos(kexact * i)

    if (deriv == 1):
        lam = abs(numodd / den)
    elif (deriv == 2):
        lam = np.sqrt(-numeven / den)

#    lam = 2.0 * np.pi * lam / N

#    plt.figure(1, figsize=(14, 12))
#    plt.plot(kexact/max(kexact), lam/kmax, linewidth=3, label='$B_i^{%d}$'%degree)
#    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
#    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=36)
#    plt.xlabel('$k/k_{max}$', fontsize=36)
#    plt.xticks(fontsize=36)
#    plt.yticks(fontsize=36)
#    plt.legend(loc='best', fontsize=36)
#    plt.tight_layout()
#    plt.show()

    return kexact, lam, kmax


def test_collocation_periodic(deriv, pts):

    degree = np.array([2, 6, 10, 14])
    degree = np.array([3, 7, 11, 15])

    plt.figure(1, figsize=(14, 12))
    for i in range(len(degree)):
        kexact, mod, kmax = collocation_periodic(deriv, degree[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, linewidth=3, label='$B_i^{%d}$'%degree[i])

    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=36)
    plt.xlabel('$k/k_{max}$', fontsize=36)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(loc='best', fontsize=36)
    plt.title('Second Derivative \n Collocation', fontsize=40)
    plt.tight_layout()
    plt.savefig('modwavenumber-second-derivative-periodic-odd-degree-collocation.eps', dpi=300)
    plt.show()


def test_compare_smooth_collocation(deriv, pts):

    lamda = [0.0, 1.0e-4/pts**2, 1.0e-5/pts**2, 1.0e-6/pts**2]
    codegree = [5]
    color = ['m', 'r', 'g', 'c', 'y', 'b']

    plt.figure(1, figsize=(14, 12))

    for i in range(len(codegree)):
        kexact, mod, kmax = collocation_smooth_periodic('original', codegree[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i], linewidth=3, label='%dth-Order B-spline'%(codegree[i]+1))

    for i in range(len(lamda)):
        kexact, mod, kmax = collocation_smooth_periodic(lamda[i], codegree[0], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i+1], \
                linewidth=3, label='%dth-Order B-spline \n(smoothed, lambda = %.1e)'%(codegree[0]+1, lamda[i]*pts**2))

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(24, 24))
    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=30)
    plt.xlabel('$k/k_{max}$', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    leg = plt.legend(bbox_to_anchor= (1.02, 0), loc='lower left', fontsize=24)
    leg.get_frame().set_alpha(0.3)
    plt.title('First Derivative \n Modified Wavenumber ', fontsize=30)
    plt.tight_layout(rect=[0, 0, 0.55, 0.65])
    plt.savefig('modwavenumber-smooth-collocation.eps', transparent=True, dpi=300)
    plt.show()


def test_compare_smooth_objective(pts):
    import matplotlib

    p_deriv = [2,3,4]
    codegree = [5]
    color = ['m', 'r', 'g', 'c', 'y', 'b']

    random = np.zeros(pts)
    for i in range(1,pts-1):
        random[i] = np.random.random(1) * 0.1

    plt.figure(1, figsize=(14, 12))
    #plt.xticks(fontsize=20)
    #plt.yticks(fontsize=20)

    matplotlib.rc('xtick', labelsize=24)
    matplotlib.rc('ytick', labelsize=24)

    for i in range(len(p_deriv)):
        x, u, l, f, lstar, fstar = collocation_smooth_periodic(0.0, codegree[0], pts, 'objective', p_deriv[i], random)
        plt.plot(l, f, color[i+1], \
                linewidth=3, label='Derivative %d penalty,  $\\lambda$ = %.2e'%(p_deriv[i], lstar))
        plt.plot(lstar, fstar, 'ro', markersize=10)

    #plt.rc('text', usetex=True)
    plt.figure(1, figsize=(24, 24))
    plt.ylabel('$V(\\lambda)$', fontsize=30)
    plt.xlabel('$\\lambda$', fontsize=30)
    plt.ylim(top=np.max(f[1:]*1.2), bottom=0.0007)
    plt.xlim(left=-0.0, right=0.1)
#    leg = plt.legend(bbox_to_anchor= (1.02, 0), loc='lower left', fontsize=24)
#    leg.get_frame().set_alpha(0.3)
    plt.legend(fontsize=24, loc='upper right')
    plt.title('Objective function plot (Degree 5 B-spline)', fontsize=30)
#    plt.tight_layout(rect=[0, 0, 0.55, 0.65])
    plt.tight_layout()
    plt.savefig('objective.eps', transparent=True, dpi=300)
    plt.show()


def test_compare_compact_collocation(deriv, pts):

    codegree = [3,  7]
    comdegree = [4, 8]
    a = [0.0, 1.0/3.0, 3.0/8.0]
    color = ['m', 'r', 'g', 'c', 'y']

    plt.figure(1, figsize=(18, 14))
    for i in range(len(codegree)):
        kexact, mod, kmax = collocation_periodic(deriv, codegree[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i], linewidth=3, label='%dth-Order B-spline'%(codegree[i]+1))

    for i in range(len(comdegree)):
        kexact, mod, kmax = compact_fd_tri(a[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i] + '--', linewidth=3, label='%dth-Order Compact FD'%comdegree[i])

    #plt.rc('text', usetex=True)
    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=30)
    plt.xlabel('$k/k_{max}$', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    leg = plt.legend(bbox_to_anchor= (1.02, 0), loc='lower left', fontsize=28)
    leg.get_frame().set_alpha(0.3)
    plt.title('First Derivative \n Modified Wavenumber ', fontsize=30)
    plt.tight_layout(rect=[0, 0, 0.55, 0.65])
    plt.savefig('modwavenumber-compact-collocation.eps', transparent=True, dpi=300)
    plt.show()


def test_article_compare_compact_collocation(deriv, pts):

    codegree = [3,  7]
    comdegree = [4, 8]
    a = [0.0, 1.0/3.0, 3.0/8.0]
    color = ['b', 'r', 'g', 'c', 'y']

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(10, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    for i in range(len(codegree)):
        kexact, mod, kmax = collocation_periodic(deriv, codegree[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i], linewidth=3, label='%dth order B-spline'%(codegree[i]+1))

    for i in range(len(comdegree)):
        kexact, mod, kmax = compact_fd_tri(a[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i] + '--', linewidth=3, label='%dth order Compact'%comdegree[i])

    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
    plt.ylabel('$\\tilde{k}/k_{\\rm{{max}}}$', fontsize=30)
    plt.xlabel('$k/k_{\\rm{max}}$', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(loc='best', fontsize=28, frameon=False)
    #plt.title('First Derivative \n Modified Wavenumber ', fontsize=30)
    plt.tight_layout()
    plt.savefig('modwavenumber-compact-collocation.eps', transparent=True, dpi=300)
    plt.show()


def test_compare_galerkin_collocation(deriv, pts):

    codegree = np.array([3, 5, 7, 9], dtype=int)
    gadegree = np.array((codegree - 1) / 2, dtype=int)
    color = ['g', 'c', 'r', 'm', 'y']

    plt.figure(1, figsize=(16, 14))
    for i in range(len(gadegree)):
        kexact, mod, kmax = galerkin_periodic(deriv, gadegree[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i] + 'o', ms=9, markeredgecolor='none', label='Galerkin %dth-Degree'%gadegree[i])

    for i in range(len(codegree)):
        kexact, mod, kmax = collocation_periodic(deriv, codegree[i], pts)
        plt.plot(kexact/max(kexact), mod/kmax, color[i], linewidth=3, label='Collocation %dth-Degree'%codegree[i])


    plt.plot(kexact/max(kexact), kexact/max(kexact),'k', linewidth=4, label='Exact')
    plt.ylabel('$\\tilde{k}/k_{max}$', fontsize=36)
    plt.xlabel('$k/k_{max}$', fontsize=36)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    leg = plt.legend(bbox_to_anchor= (1.04, 0), loc='lower left', fontsize=25)
    leg.get_frame().set_alpha(0.3)
    plt.title('First Derivative', fontsize=40)
    plt.tight_layout(rect=[0, 0, 0.65, 0.65])
    plt.savefig('modwavenumber-galerkin-collocation.pdf', transparent=True, dpi=300)
    plt.show()


def test_compare_all(deriv, pts):
    import scipy.stats.mstats as stats

#    codegree = np.array([5], dtype=int)
#    gadegree = np.array((codegree - 1) / 2, dtype=int)
#    comdegree = [6]
#    a = [1./3.]

    codegree = np.array([5, 7, 9, 11], dtype=int)
    gadegree = np.array((codegree - 1) / 2, dtype=int)
#    comdegree = [4]
#    a = [1./4.]

    color = ['g', 'c', 'r', 'm', 'y']

    plt.figure(1, figsize=(18, 18))
    for i in range(len(gadegree)):
        kexact, mod, kmax = galerkin_periodic(deriv, gadegree[i], pts)
        sol = kexact / max(kexact)
        num = mod / kmax
        err = abs(sol - num)

        if deriv == 1:
            x = np.power(10.0, -14)
            y = 5

        elif deriv == 2:
            if codegree[i] == 5:
                x = np.power(10.0, -3)
                y = 50
            elif codegree[i] == 7:
                x = np.power(10.0, -6)
                y = 50
            elif codegree[i] == 9:
                x = np.power(10.0, -9)
                y = 40
            else:
                x = np.power(10.0, -9)
                y = 100

        idx = np.where(err > x)[0][0]
        err = np.nan_to_num(err)
        order, _, _, _ = stats.theilslopes(np.log(err[idx:(idx+y)]), np.log(sol[idx:(idx+y)]))
        plt.loglog(sol, err, color[i] + 'o', ms=9, markeredgecolor='none', \
                label='Galerkin %d-bands,   Order = %0.3f'%(codegree[i], order))

        err = np.nan_to_num(sol - num)
        err = abs(err) /(2.0 * np.pi / pts)
        idx1 = min(np.where(err >= 10.0/100.0)[0])
        idx2 = min(np.where(err >= 1.0/100.0)[0])
        idx3 = min(np.where(err >= 0.1/100.0)[0])
        print(codegree[i], '&', '%0.2f'%sol[idx1-1], '&', '%0.2f'%sol[idx2-1], '&', '%0.2f'%sol[idx3-1], '\\\\')


    for i in range(len(codegree)):
        kexact, mod, kmax = collocation_periodic(deriv, codegree[i], pts)
        sol = kexact / max(kexact)
        num = mod / kmax
        err = abs(sol - num)
        plt.loglog(sol/kmax, err/kmax, color[i], linewidth=3, label='Collocation %d-bands'%codegree[i])

#        err = np.nan_to_num(sol - num)
#        print('Collocation')
#        idx = min(np.where(err >= 10.0/100.0)[0])
#        print(codegree[i], sol[idx-1])
#        idx = min(np.where(err >= 1.0/100.0)[0])
#        print(codegree[i], sol[idx-1])
#        idx = min(np.where(err >= 0.1/100.0)[0])
#        print(codegree[i], sol[idx-1])


#    for i in range(len(comdegree)):
#        kexact, mod, kmax = compact_fd_tri(a[i], pts)
#        sol = kexact / max(kexact)
#        num = mod / kmax
#        err = sol - num
#        plt.loglog(sol, err, color[i] + '--', linewidth=3, label='%dth-Order Compact'%comdegree[i])

    # plt.plot(sol, sol,'k', linewidth=4, label='Exact')

    plt.ylabel('Error', fontsize=36)
    plt.xlabel('$k/k_{max}$', fontsize=36)
    plt.xlim(1e-2, 1e-0)
    plt.ylim(1e-12, 1e-0)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    leg = plt.legend(bbox_to_anchor= (1.04, 0), loc='lower left', fontsize=25)
    leg.get_frame().set_alpha(0.3)
    plt.tight_layout(rect=[0, 0, 0.55, 0.6])
    if deriv == 1:
        plt.title('First Derivative', fontsize=40)
        plt.savefig('modwavenumber-first-derivative-all-error.eps', transparent=True, dpi=300)
    elif deriv == 2:
        plt.title('Second Derivative', fontsize=40)
        plt.savefig('modwavenumber-second-derivative-all-error.eps', transparent=True, dpi=300)
    plt.show()


def main():
    # collocation_smooth_periodic(0.0, 5, 32, 'objective', 2)
    # galerkin_periodic(1, 5, 200)
    # compact_fd_tri(1.0/3.0, 200)
    # collocation_periodic(1, 5, 200)
    # test_collocation_periodic(2, 200)
    #test_compare_compact_collocation(1, 200)
    test_article_compare_compact_collocation(1, 200)
    # test_compare_galerkin_collocation(1, 150)
    # test_compare_smooth_collocation(1, 50)
    # test_compare_smooth_objective(32)
    # test_compare_all(2, 200)


if __name__ == "__main__":
    main()
