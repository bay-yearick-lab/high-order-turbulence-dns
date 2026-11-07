import numpy as np
import matplotlib.pyplot as plt


def test_channel_flow_ic():

    mink=0.0; maxk=2.0
    pts = 200
    u0 = np.zeros(pts)
    x2 = np.linspace(mink, maxk, pts)

    # plot
    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(9, 7))
    plt.rcParams["font.family"] = "Times New Roman"

    u = 1.5 * (1.0 - (x2-1.0)**2)
    T = 1.0 + (0.4) * 0.72 * 1.5**2 / 3.0 * 1.5 * (1.0 - (x2-1.0)**4)

    plt.plot(u, x2, 'b--', linewidth=3, label='$u_1(t=0)$')
    plt.fill_between(u, u0, x2, color='b', alpha=0.2)

    plt.plot(T, x2, 'r--', linewidth=3, label='$T(t=0, T_w = 1)$')
    plt.fill_between(T, u0, x2, color='r', alpha=0.2)

    plt.ylim(0, 2.0)
    plt.xlim(-0.1, 1.6)
    plt.xticks([0.0, 0.5, 1.0, 1.5], fontsize=30)
    plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0], fontsize=30)
    plt.ylabel('$x_2$', labelpad=8, fontsize=35)
    plt.legend(loc=6, fontsize=35, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('channel_ic_viscosity.pdf', transparent=True, dpi=300)

    plt.show()
    return

def plot_bspline_piecewise():

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(9, 8))
    for j in range(-3, 3):
        u = np.linspace(j, j+1, 100)
        y = np.zeros(len(u))
        for i in range(len(u)):
            x = u[i]
            if -3<x<=-2:
                y[i] = 1/120 * (3+x)**5
            elif -2<x<=-1:
                y[i] = 1/120 * (51-75*x-210*x**2-150*x**3-45*x**4-5*x**5)
            elif -1<x<=0:
                y[i] = 1/60 * (33-30*x**2 + 15*x**4 +5*x**5)
            elif 0<x<=1:
                y[i] = 1/60 * (33-30*x**2 + 15*x**4 -5*x**5)
            elif 1<x<=2:
                y[i] = 1/120 * (51+75*x-210*x**2+150*x**3-45*x**4+5*x**5)
            elif 2<x<=3:
                y[i] = -1/120 * (-3+x)**5

        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.title("B-spline degree 5", fontsize=40)
        plt.plot(u, y, linewidth=7)
        plt.savefig("bspline-degree5-piecewise.eps", transparent=True, dpi=300)
    plt.show()

def test_analytic_bspline_svd(degree, pts):
    import pygsvd
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]

    _ ,_ ,_ ,_ , fullP_1D = \
    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    print(fullP_1D)
    quit()

    row, col = fullB_1D.shape
    periodic_row = np.zeros(row)
    periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

    periodic_row = np.roll(periodic_row, -2)

    for i in range(row):
        fullB_1D[i, :] = np.roll(periodic_row, i)

    AAT = fullB_1D @ fullB_1D.T
    band = AAT[0, :]
    eigenvalue = np.zeros(len(band), dtype=complex)
    w = np.zeros(len(band), dtype=complex)

    for i in range(len(band)):
        w[i] = np.exp(2.0 * np.pi * i /len(band) * 1j)

    for j in range(len(band)):
        for i in range(len(band)):
            eigenvalue[j] = eigenvalue[j] + band[i] * w[j]**i

    eigenvalue = np.real(np.sqrt(eigenvalue))
    idx = eigenvalue.argsort()[::-1]
    eigenvalue = eigenvalue[idx]

    eigenvector = np.zeros((len(band), len(band)), dtype=complex)
    left_eigenvector = np.zeros((len(band), len(band)))
    for j in range(len(band)):
        for i in range(len(band)):
            eigenvector[j, i] = w[j]**i / np.sqrt(len(band))

    left_eigenvector = np.real(eigenvector[:, idx])

    U, S, VT = la.svd(fullB_1D)

    # check singular values
    print(np.allclose(eigenvalue, S))

    x = np.linspace(0, len(band), len(band))

    #V = left_eigenvector
    V = U

    plt.rc('text', usetex=True)
    fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(20, 10))
    for i in range(10):
        plt.subplot(2, 5, i+1)
        plt.plot(x, np.real(V[:,i]))
        plt.ylim(top=0.35, bottom=-0.35)
        plt.xticks(np.arange(0, pts, step=50), fontsize=15)
        plt.yticks(np.arange(-0.2, 0.22, step=0.2), fontsize=15)

        plt.title('$i = {%d}, ~\sigma_{%d} = %0.4f$'%(i+1, i+1, S[i].real), \
                fontsize=20, bbox=dict(facecolor='none', edgecolor='black', pad=5.0), position=(0.58, 1.02))
    #plt.savefig('allv.eps', transparent=True, dpi=300)

    fig1, ax1 = plt.subplots(nrows=2, ncols=5, figsize=(20, 10))
    for i in range(10):
        plt.subplot(2, 5, i+1)
        plt.plot(x, np.real(V[:,-(i+1)]), \
                label='$i = %d$, $\sigma = %.4f$'%(V.shape[0]-i, S[-(i+1)].real))
        plt.ylim(top=0.35, bottom=-0.35)
        plt.xticks(np.arange(0, pts, step=50), fontsize=15)
        plt.yticks(np.arange(-0.2, 0.22, step=0.2), fontsize=15)

        plt.title('$i = {%d}, ~\sigma_{%d} = %0.4f$'%(V.shape[0]-i, V.shape[0]-i, S[-(i+1)].real), \
                fontsize=20, bbox=dict(facecolor='none', edgecolor='black', pad=5.0), position=(0.56, 1.02))

    #plt.savefig('allvbehind.eps', transparent=True, dpi=300)
    plt.show()

def test_ideal_transfer():

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    pts = 200
    kexact = np.linspace(0.0, np.pi, pts)
    ideal = np.zeros(pts)

    for i in range(pts):
        if kexact[i] <  0.75 * np.pi:
            ideal[i] = 1.0
        else:
            ideal[i] = 0.0

    plt.plot(kexact, ideal, 'k', linewidth=3)
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.tight_layout()
    plt.savefig('transfer-ideal.eps', transparent=True, dpi=300)
    plt.show()
    return


def test_boundary_fd_transfer(pts):

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi

    kexact = np.linspace(0.0, np.pi, pts)
    compactTF6 = np.linspace(mink, maxk, pts, dtype=complex)
    compactTF8 = np.linspace(mink, maxk, pts, dtype=complex)

    # compact schemes
    alpha = 0.49
    a62 = np.array([1./64.+31*alpha/32., 29./32.+3.*alpha/16., 15./64.+17.*alpha/32, -5./16.+5*alpha/8., \
          15./64.-15*alpha/32., -3./32+3*alpha/16, 1./64.-alpha/32.])
    a82 = np.array([1./256 + 127*alpha/128., 31/32 + alpha/16., 7./64 + 25*alpha/32., -7/32. + 7*alpha/16,\
            7*(5-10.*alpha)/128, -7/32. + 7.*alpha/16, 7/64.-7*alpha/32, -1/32.+alpha/16., 1/256.-alpha/128])
    compactlabel=['6th order boundary compact ($\\alpha=0.49$)',\
                  '8th order boundary compact ($\\alpha=0.49$)']
    compactcolor=['b--', 'r--']

    den = 1.0 + 2.0 * alpha * np.cos(kexact)
    for i in range(pts):
        num = 0.0
        for n in range(len(a62)):
            num = num  + a62[n] * (np.cos((n-1) * kexact[i]) + 1j*np.sin((n-1) * kexact[i]))
        compactTF6[i] = num / den[i]

    den = 1.0 + 2.0 * alpha * np.cos(kexact)
    for i in range(pts):
        num = 0.0
        for n in range(len(a82)):
            num = num  + a82[n] * (np.cos((n-1) * kexact[i]) + 1j*np.sin((n-1) * kexact[i]))
        compactTF8[i] = num / den[i]

    drpTF2 = np.linspace(mink, maxk, pts, dtype=complex)

    # drp schemes
    a27 = np.array([-0.057717512738, 0.199278374994, -0.292668277650, \
            0.244537361546, -0.134605018019, 0.056184263460, -0.015009191593])
    drplabel=['2nd order boundary optimized finite-difference']
    drpcolor=['g--']

    for i in range(pts):
        num = 0.0
        for n in range(len(a62)):
            num = num  + a27[n] * (np.cos((n-1) * kexact[i]) + 1j*np.sin((n-1) * kexact[i]))
        drpTF2[i] = 1.0 - num

    # plot real
    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"
    plt.plot(kexact, np.real(compactTF6), compactcolor[0], linewidth=3, label=compactlabel[0])
    plt.plot(kexact, np.real(compactTF8), compactcolor[1], linewidth=3, label=compactlabel[1])
    plt.plot(kexact, np.real(drpTF2), drpcolor[0], linewidth=3, label=drplabel[0])
    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(-0.1, 1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (real)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.legend(fontsize=28, framealpha=0.0)
    #plt.savefig('transfer-real-boundary-compact-bspline.pdf', transparent=True, dpi=300)
    plt.tight_layout()

    # plot imag
    plt.rc('text', usetex=True)
    plt.figure(2, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"
    plt.plot(kexact, np.imag(compactTF6), compactcolor[0], linewidth=3, label=compactlabel[0])
    plt.plot(kexact, np.imag(compactTF8), compactcolor[1], linewidth=3, label=compactlabel[1])
    plt.plot(kexact, np.imag(drpTF2), drpcolor[0], linewidth=3, label=drplabel[0])
    plt.plot(kexact, np.zeros(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(-0.1, 1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (imaginary)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.legend(fontsize=28, framealpha=0.0)
    #plt.savefig('transfer-imag-boundary-compact-bspline.pdf', transparent=True, dpi=300)

    plt.show()
    return


def output_derivative_order(degree):
    import csv

    rhoerror = []
    uerror = []
    verror = []
    werror = []
    Terror = []
    pts = np.array([40, 60, 80, 100, 120, 140, 160, 180])

    for i in range(len(pts)):
        print(pts[i])
        a, b, c, d, e = test_derivative_order(degree, pts[i])

        rhoerror.append(a)
        uerror.append(b)
        verror.append(c)
        werror.append(d)
        Terror.append(e)

        temp = zip(pts, rhoerror, uerror, verror, werror, Terror)

        with open('derivative-degree%d.csv'% degree, "w") as output:
            writer = csv.writer(output, lineterminator='\n', delimiter=',')
            for val in temp:
                writer.writerow(val)
    return


def test_derivative_order(degree, pts):
    import sympy as sp
    import numpy.linalg as la
    import scipy.linalg as sla
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 0.3
    # number of derivatives
    deriv = 2

    # B-spline
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)
    B = ders[:, 0, :]
    L = ders[:, 2, :]
    LU = sla.lu_factor(B)

    # sympy rho
    xx = sp.symbols('xx')
    rho = sp.exp(sp.sin(6.0 * np.pi * xx)) * sp.exp(sp.sin(8.0 * np.pi * xx)) * sp.exp(sp.sin(6.0 * np.pi * xx)) + 2.0
    rho_xx = sp.lambdify((xx), rho.diff(xx, xx), 'numpy')
    rho = sp.lambdify((xx), rho, 'numpy')

    xx = sp.symbols('xx')
    u = sp.exp(sp.sin(8.0 * np.pi * xx)) * sp.exp(sp.cos(6.0 * np.pi * xx)) * sp.exp(sp.cos(8.0 * np.pi * xx))
    u_xx = sp.lambdify((xx), u.diff(xx, xx), 'numpy')
    u = sp.lambdify((xx), u, 'numpy')

    xx = sp.symbols('xx')
    v = sp.exp(sp.cos(6.0 * np.pi * xx)) * sp.exp(sp.sin(8.0 * np.pi * xx)) * sp.exp(sp.cos(6.0 * np.pi * xx))
    v_xx = sp.lambdify((xx), v.diff(xx, xx), 'numpy')
    v = sp.lambdify((xx), v, 'numpy')

    xx = sp.symbols('xx')
    w = sp.exp(sp.cos(8.0 * np.pi * xx)) * sp.exp(sp.cos(6.0 * np.pi * xx)) * sp.exp(sp.sin(6.0 * np.pi * xx))
    w_xx = sp.lambdify((xx), w.diff(xx, xx), 'numpy')
    w = sp.lambdify((xx), w, 'numpy')

    xx = sp.symbols('xx')
    T = sp.exp(sp.cos(6.0 * np.pi * xx)) * sp.exp(sp.cos(8.0 * np.pi * xx)) * sp.exp(sp.cos(6.0 * np.pi * xx)) + 2.0
    T_xx = sp.lambdify((xx), T.diff(xx, xx), 'numpy')
    T = sp.lambdify((xx), T, 'numpy')

    # rho
    alpha = sla.lu_solve(LU, rho(fullu))
    num = L @ alpha
    rhoerror = la.norm(rho_xx(fullu) - num, np.inf)

    # u
    alpha = sla.lu_solve(LU, u(fullu))
    num = L @ alpha
    uerror = la.norm(u_xx(fullu) - num, np.inf)

    # v
    alpha = sla.lu_solve(LU, v(fullu))
    num = L @ alpha
    verror = la.norm(v_xx(fullu) - num, np.inf)

    # W
    alpha = sla.lu_solve(LU, w(fullu))
    num = L @ alpha
    werror = la.norm(w_xx(fullu) - num, np.inf)

    # T
    alpha = sla.lu_solve(LU, T(fullu))
    num = L @ alpha
    Terror = la.norm(T_xx(fullu) - num, np.inf)

    return rhoerror, uerror, verror, werror, Terror


def test_eig_nonperiodic_drp(degree, pts):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    kexactB = np.linspace(0.0, np.pi, pts-4)
    kexact = np.linspace(0.0, np.pi, pts)

    Erow = np.zeros(pts)

    # 9 POINT STENCIL
    Erow.fill(0.0)
    # 4th order compact (DRP)
    a = 0.243527493120
    b = -0.204788880640
    c = 0.120007591680
    d = -0.045211119360
    e = 0.008228661760

    Erow.fill(0.0)

    Erow[-4] = e; Erow[-3] = d; Erow[-2] = c; Erow[-1] = b; Erow[0] = a
    Erow[4] = e; Erow[3] = d; Erow[2] = c; Erow[1] = b

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    for j in range(1):
        if j==0:
            a42 = np.array([-0.057717512738,0.199278374994,-0.292668277650,0.244537361546,\
                    -0.134605018019,0.056184263460,-0.015009191593, 0.0, 0.0, 0.0, 0.0])
            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0])
            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090])
        elif j==1:
            a42 = np.array([-0.085777408970,0.277628171524,-0.356848072173,0.223119093072,\
                    -0.057347064865,-0.000747264596,-0.000027453993, 0.0, 0.0, 0.0, 0.0])
            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0])
            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090])

        B4up = np.zeros((3, pts))
        for i in range(len(a43)):
            B4up[0, i] = a42[i]
            B4up[1, i] = a43[i]
            B4up[2, i] = a44[i]

        B4down = np.fliplr(np.flipud(B4up))
        explicit4 = np.zeros((pts, pts))

        for i in range(pts):
            explicit4[i, :] = np.roll(Erow, i)

        explicit4[0, :]  = 0.0; explicit4[0, 0]   = 1.0
        explicit4[-1, :] = 0.0; explicit4[-1, -1] = 1.0
        explicit4[1:4, :] = B4up
        explicit4[-4:-1, :] = B4down

        eig4 = la.eigvals(np.eye(pts)-explicit4)

        print('Scheme %i, 4th order (9 point), minreal=%.6f, maxreal=%.6f, minimag=%.6f, maximag=%.6f'
            %(j, np.min(np.real(eig4)), np.max(np.real(eig4)), np.min(np.imag(eig4)), np.max(np.imag(eig4))))
        print('minMag=%.6f, maximag=%.6f'%(np.min(np.absolute(eig4)), np.max(np.absolute(eig4))))

        # plot real
        plt.plot(kexact, np.sort(np.real(eig4))[::-1], 'r--', linewidth=3, \
        label='Optimized 9pt-\\textrm{\MakeUppercase{\\romannumeral %d}}'%(j+1))

#    # 11 POINT STENCIL
#    Erow.fill(0.0)
#    # 2nd order compact (DRP)
#    a = 0.215044884112
#    b = -0.187772883589
#    c = 0.123755948787
#    d = -0.059227575576
#    e = 0.018721609157
#    f = -0.002999540835
#
#    Erow.fill(0.0)
#
#    Erow[-5] = f; Erow[-4] = e; Erow[-3] = d; Erow[-2] = c; Erow[-1] = b; Erow[0] = a
#    Erow[5] = f; Erow[4] = e; Erow[3] = d; Erow[2] = c; Erow[1] = b
#
#    for j in range(2):
#        if j==0:
#            a42 = np.array([-0.057717512738,0.199278374994,-0.292668277650,0.244537361546,\
#                    -0.134605018019,0.056184263460,-0.015009191593, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
#                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
#            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090, 0.0, 0.0])
#            a45 = np.array([0.008391235145,-0.047402506444,0.121438547725,-0.200063042812,0.240069047836,\
#                    -0.207269200140,0.122263107844,-0.047121062819,0.009014891495,0.001855812216,-0.001176830044, 0.0, 0.0])
#        elif j==1:
#            a42 = np.array([-0.085777408970,0.277628171524,-0.356848072173,0.223119093072,\
#                    -0.057347064865,-0.000747264596,-0.000027453993, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
#                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
#            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090, 0.0, 0.0])
#            a45 = np.array([0.008391235145,-0.047402506444,0.121438547725,-0.200063042812,0.240069047836,\
#                    -0.207269200140,0.122263107844,-0.047121062819,0.009014891495,0.001855812216,-0.001176830044, 0.0, 0.0])
#
#        B4up = np.zeros((4, pts))
#        for i in range(len(a43)):
#            B4up[0, i] = a42[i]
#            B4up[1, i] = a43[i]
#            B4up[2, i] = a44[i]
#            B4up[3, i] = a45[i]
#
#        B4down = np.fliplr(np.flipud(B4up))
#        explicit4 = np.zeros((pts, pts))
#
#        for i in range(pts):
#            explicit4[i, :] = np.roll(Erow, i)
#
#        explicit4[0, :]  = 0.0; explicit4[0, 0]   = 1.0
#        explicit4[-1, :] = 0.0; explicit4[-1, -1] = 1.0
#        explicit4[1:5, :] = B4up
#        explicit4[-5:-1, :] = B4down
#
#        eig4 = la.eigvals(np.eye(pts)-explicit4)
#
#        print('Scheme %i, 2nd order (11 point), minreal=%.6f, maxreal=%.6f, minimag=%.6f, maximag=%.6f'\
#        %(j, np.min(np.real(eig4)), np.max(np.real(eig4)), np.min(np.imag(eig4)), np.max(np.imag(eig4))))
#        print('minMag=%.6f, maximag=%.6f'%(np.min(np.absolute(eig4)), np.max(np.absolute(eig4))))
#
#        # plot real
#        plt.plot(kexact, np.sort(np.real(eig4))[::-1], linestyle='--', linewidth=3, \
#        label='Optimized 11pt, \\textrm{\MakeUppercase{\\romannumeral %d}}'%(j+1))
#
#    # 13 POINT STENCIL
#    Erow.fill(0.0)
#    # 4th order compact (DRP)
#    a = 0.215044884112
#    b = -0.187772883589
#    c = 0.123755948787
#    d = -0.059227575576
#    e = 0.018721609157
#    f = -0.002999540835
#    g = -0.002999540835
#
#
#    Erow.fill(0.0)
#
#    Erow[-6] = g; Erow[-5] = f; Erow[-4] = e; Erow[-3] = d; Erow[-2] = c; Erow[-1] = b; Erow[0] = a
#    Erow[6] = g; Erow[5] = f; Erow[4] = e; Erow[3] = d; Erow[2] = c; Erow[1] = b
#
#    for j in range(2):
#        if j==0:
#            a42 = np.array([-0.057717512738,0.199278374994,-0.292668277650,0.244537361546,\
#                    -0.134605018019,0.056184263460,-0.015009191593, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
#                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
#            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090, 0.0, 0.0])
#            a45 = np.array([0.008391235145,-0.047402506444,0.121438547725,-0.200063042812,0.240069047836,\
#                    -0.207269200140,0.122263107844,-0.047121062819,0.009014891495,0.001855812216,-0.001176830044, 0.0, 0.0])
#            a46 = np.array([-0.002999540835,0.018721609157,-0.059227575576,0.123755948787,-0.187772883589,0.215044884112,\
#                    -0.187772883589,0.123755948787,-0.059227575576,0.018721609157,-0.002999540835, 0.0, 0.0])
#
#        elif j==1:
#            a42 = np.array([-0.085777408970,0.277628171524,-0.356848072173,0.223119093072,\
#                    -0.057347064865,-0.000747264596,-0.000027453993, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
#                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
#            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090, 0.0, 0.0])
#            a45 = np.array([0.008391235145,-0.047402506444,0.121438547725,-0.200063042812,0.240069047836,\
#                    -0.207269200140,0.122263107844,-0.047121062819,0.009014891495,0.001855812216,-0.001176830044, 0.0, 0.0])
#            a46 = np.array([-0.002999540835,0.018721609157,-0.059227575576,0.123755948787,-0.187772883589,0.215044884112,\
#                    -0.187772883589,0.123755948787,-0.059227575576,0.018721609157,-0.002999540835, 0.0, 0.0])
#
#
#        B4up = np.zeros((5, pts))
#        for i in range(len(a43)):
#            B4up[0, i] = a42[i]
#            B4up[1, i] = a43[i]
#            B4up[2, i] = a44[i]
#            B4up[3, i] = a45[i]
#            B4up[4, i] = a46[i]
#
#        B4down = np.fliplr(np.flipud(B4up))
#        explicit4 = np.zeros((pts, pts))
#
#        for i in range(pts):
#            explicit4[i, :] = np.roll(Erow, i)
#
#        explicit4[0, :]  = 0.0; explicit4[0, 0]   = 1.0
#        explicit4[-1, :] = 0.0; explicit4[-1, -1] = 1.0
#        explicit4[1:6, :] = B4up
#        explicit4[-6:-1, :] = B4down
#
#        eig4 = la.eigvals(np.eye(pts)-explicit4)
#
#        print('Scheme %i, 4th order (13 point), minreal=%.6f, maxreal=%.6f, minimag=%.6f, maximag=%.6f'\
#        %(j, np.min(np.real(eig4)), np.max(np.real(eig4)), np.min(np.imag(eig4)), np.max(np.imag(eig4))))
#        print('minMag=%.6f, maximag=%.6f'%(np.min(np.absolute(eig4)), np.max(np.absolute(eig4))))
#
#        # plot real
#        plt.plot(kexact, np.sort(np.real(eig4))[::-1], linestyle='--', linewidth=3, \
#        label='Optimized 13pt, \\textrm{\MakeUppercase{\\romannumeral %d}}'%(j+1))

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    kexactB = np.linspace(0.0, np.pi, pts-4)
    kexact = np.linspace(0.0, np.pi, pts)

    # B-spline
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)
    B = ders[:, 0, :]
    I = np.eye(B.shape[0])
#    _ ,_ ,_ ,_ , P = \
#    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, B, 2, 'nonperiodic', matrix_type='full')

    A = (B @ la.inv(B.T @ B + 1e-8*I) @ B.T)[1:-1, 1:-1]
    eig = la.eigvals(A)

    eig = np.sort(np.real(eig))[::-1]
    eig = (eig - np.min(eig))/np.ptp(eig)
    eig = eig[:-2]
    eig = (eig - np.min(eig))/np.ptp(eig)

    plt.plot(kexactB, (np.real(eig)), 'b--',  linewidth=3,  \
            label='8th-order B-spline ($B^{(6)}, \lambda^2 = 10^{-6}$)')
    plt.plot(kexact, np.ones(len(kexact)), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=np.min(np.real(eig4))-0.1, top=1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('tf-real-optimized-8-bspline.eps', transparent=True, dpi=300)

    plt.show()
    return


def test_eig_nonperiodic_drp_7point_antisym(degree, pts):
    import numpy.linalg as la
    import scipy.linalg as sla

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    kexactB = np.linspace(0.0, np.pi, pts-4)
    kexact = np.linspace(0.0, np.pi, pts)

    Erow = np.zeros(pts)

    for k in range(5):
        # 6th order compact (DRP)
        if k==0:
            a = -0.050000000000
            b = 0.246093750000
            c = -0.450000000000
            d = 0.328125000000
            e = -0.098437500000
            f = 0.028125000000
            g = -0.003906250000
        elif k==1:
            a = 0.031250000000
            b = 0.000000000000
            c = -0.156250000000
            d = 0.312500000000
            e = -0.312500000000
            f = 0.156250000000
            g = -0.031250000000
        elif k==2:
            a = 0.000000000000
            b = 0.031250000000
            c = -0.156250000000
            d = 0.312500000000
            e = -0.312500000000
            f = 0.156250000000
            g = -0.031250000000
        elif k==3:
            a = -0.027343750000
            b = 0.150000000000
            c = -0.328125000000
            d = 0.350000000000
            e = -0.164062500000
            f = 0.021875000000
            g = -0.002343750000
        elif k==4:
            a = -0.018750000000
            b = 0.109375000000
            c = -0.262500000000
            d = 0.328125000000
            e = -0.218750000000
            f = 0.065625000000
            g = -0.003125000000

        Erow.fill(0.0)

        Erow[-3] = a; Erow[-2] = b; Erow[-1] = c; Erow[0] = d
        Erow[3] = e; Erow[2] = f; Erow[1] = g

        for j in range(2):
            if j==0:
                a42 = np.array([-0.057717512738,0.199278374994,-0.292668277650,0.244537361546,\
                        -0.134605018019,0.056184263460,-0.015009191593, 0.0, 0.0, 0.0, 0.0])
                a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
                        0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0])

            elif j==1:
                a42 = np.array([-0.085777408970,0.277628171524,-0.356848072173,0.223119093072,\
                        -0.057347064865,-0.000747264596,-0.000027453993, 0.0, 0.0, 0.0, 0.0])
                a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
                        0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0])

            B4up = np.zeros((2, pts))
            for i in range(len(a43)):
                B4up[0, i] = a42[i]
                B4up[1, i] = a43[i]

            B4down = np.fliplr(np.flipud(B4up))
            explicit4 = np.zeros((pts, pts))

            for i in range(pts):
                explicit4[i, :] = np.roll(Erow, i)

            explicit4[0, :]  = 0.0; explicit4[0, 0]   = 1.0
            explicit4[-1, :] = 0.0; explicit4[-1, -1] = 1.0
            explicit4[1:3, :] = B4up
            explicit4[-3:-1, :] = B4down

            eig4 = la.eigvals(np.eye(pts)-explicit4)

            print('Scheme %i, 6th order optimized finite-difference, minreal=%.6f, maxreal=%.6f, minimag=%.6f, maximag=%.6f'\
            %(k, np.min(np.real(eig4)), np.max(np.real(eig4)), np.min(np.imag(eig4)), np.max(np.imag(eig4))))

            plt.rc('text', usetex=True)
            plt.figure(1, figsize=(16, 8))
            plt.rcParams["font.family"] = "Times New Roman"

            # plot real
            plt.plot(kexact, np.sort(np.real(eig4))[::-1], linestyle='--', linewidth=3, \
            label='4th order optimized finite-difference')

    plt.plot(kexact, np.ones(len(kexact)), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=np.min(np.real(eig4))-0.1, top=1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (real)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    #plt.savefig('tf-real-compact-8-bspline.eps', transparent=True, dpi=300)

    plt.show()
    return


def test_eig_nonperiodic_drp_9point_centered(degree, pts):
    import numpy.linalg as la
    import scipy.linalg as sla

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    kexactB = np.linspace(0.0, np.pi, pts-4)
    kexact = np.linspace(0.0, np.pi, pts)

    Erow = np.zeros(pts)
    # 4th order compact (DRP)
    a = 0.243527493120
    b = -0.204788880640
    c = 0.120007591680
    d = -0.045211119360
    e = 0.008228661760

    Erow.fill(0.0)

    Erow[-4] = e; Erow[-3] = d; Erow[-2] = c; Erow[-1] = b; Erow[0] = a
    Erow[4] = e; Erow[3] = d; Erow[2] = c; Erow[1] = b

    for j in range(2):
        if j==0:
            a42 = np.array([-0.057717512738,0.199278374994,-0.292668277650,0.244537361546,\
                    -0.134605018019,0.056184263460,-0.015009191593, 0.0, 0.0, 0.0, 0.0])
            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0])
            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090])

        elif j==1:
            a42 = np.array([-0.085777408970,0.277628171524,-0.356848072173,0.223119093072,\
                    -0.057347064865,-0.000747264596,-0.000027453993, 0.0, 0.0, 0.0, 0.0])
            a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
                    0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0])
            a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
            0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090])

        B4up = np.zeros((3, pts))
        for i in range(len(a44)):
            B4up[0, i] = a42[i]
            B4up[1, i] = a43[i]
            B4up[2, i] = a44[i]

        B4down = np.fliplr(np.flipud(B4up))
        explicit4 = np.zeros((pts, pts))

        for i in range(pts):
            explicit4[i, :] = np.roll(Erow, i)

        explicit4[0, :]  = 0.0; explicit4[0, 0]   = 1.0
        explicit4[-1, :] = 0.0; explicit4[-1, -1] = 1.0
        explicit4[1:4, :] = B4up
        explicit4[-4:-1, :] = B4down

        eig4 = la.eigvals(np.eye(pts)-explicit4)

        print('4th order optimized finite-difference, maxreal=%.15f, minreal=%.15f, maximag=%.15f, minimag=%.15f'\
        %(np.max(np.real(eig4)), np.min(np.real(eig4)), np.max(np.imag(eig4)), np.min(np.imag(eig4))))
    return


def test_eig_nonperiodic_transfer(degree, pts):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    kexactB = np.linspace(0.0, np.pi, pts-4)
    kexact = np.linspace(0.0, np.pi, pts)

    # B-spline
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)
    B = ders[:, 0, :]
    I = np.eye(B.shape[0])
#    _ ,_ ,_ ,_ , P = \
#    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, B, 2, 'nonperiodic', matrix_type='full')

    A = (B @ la.inv(B.T @ B + 1e-8*I) @ B.T)[1:-1, 1:-1]
    eig = la.eigvals(A)

    eig = np.sort(np.real(eig))[::-1]
    eig = (eig - np.min(eig))/np.ptp(eig)
    eig = eig[:-2]
    eig = (eig - np.min(eig))/np.ptp(eig)
    print('8th order B-spline, max=%.15f, min=%.15f'%(np.max(np.real(eig)), np.min(np.real(eig))))

    # 8th order compact (high order one-sided)
    alpha=0.49
    a = 1/128.0 * (93.0 + 70 * alpha)
    b = 1.0/16.0 * (7.0 + 18.0 * alpha)
    c = 1.0/32.0 * (-7.0 + 14.0 *alpha)
    d = 1.0/16.0 * (1.0 - 2.0 * alpha)
    e = 1.0/128.0 * (-1.0 + 2.0 *alpha)

    Erow = np.zeros(pts)
    Erow[-4] = e/2; Erow[-3] = d/2; Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a
    Erow[4] = e/2; Erow[3] = d/2; Erow[2] = c/2; Erow[1] = b/2
    a82 = np.array([1./256 + 127*alpha/128., 31/32 + alpha/16., 7./64 + 25*alpha/32., -7/32. + 7*alpha/16.,\
            7*(5-10.*alpha)/128., -7/32. + 7.*alpha/16, 7/64.-7*alpha/32., -1/32.+alpha/16., 1/256.-alpha/128.])
    a83 = np.array([-1./256 + alpha/128., 1/32 + 15*alpha/16., 57/64. + 7.*alpha/32. ,7./32 + 9*alpha/16., \
            -7*(5-10.*alpha)/128., 7/32. - 7.*alpha/16, -7/64.+7*alpha/32, 1/32.-alpha/16., -1/256.+alpha/128.])
    a84 = np.array([1./256 - alpha/128., -1/32 + alpha/16., 7/64. + 25.*alpha/32. ,25./32 + 7*alpha/16., \
            35/128.+29.*alpha/64, -7/32. + 7.*alpha/16, 7/64.-7*alpha/32., -1/32.+alpha/16., 1/256.-alpha/128.])

    B8up = np.zeros((3, pts))
    for i in range(len(a82)):
        B8up[0, i] = a82[i]
        B8up[1, i] = a83[i]
        B8up[2, i] = a84[i]

    B8down = np.fliplr(np.flipud(B8up))
    explicit8 = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit8 = np.zeros((pts, pts))

    for i in range(pts):
        explicit8[i, :] = np.roll(Erow, i)
        implicit8[i, :] = np.roll(Irow, i)

    explicit8[0, :]  = 0.0; explicit8[0, 0]   = 1.0
    explicit8[-1, :] = 0.0; explicit8[-1, -1] = 1.0
    implicit8[0, :]  = 0.0; implicit8[0, 0]   = 1.0
    implicit8[-1, :] = 0.0; implicit8[-1, -1] = 1.0

    explicit8[:3, :] = B8up
    explicit8[-3:, :] = B8down

    A8 = la.solve(implicit8, explicit8)
    eig8 = la.eigvals(A8)
    imag = np.imag(eig8)[::-1]
    pos_imag_eig8 = imag[imag>=0]
    kexact_pos =kexact[np.where(imag>=0)]
    neg_imag_eig8 = imag[imag<=0]
    kexact_neg =kexact[np.where(imag<=0)]

    print('8th order compact, max=%.15f, min=%.15f'%(np.max(np.real(eig8)), np.min(np.real(eig8))))

    # 8th order compact (low order centered)
    alpha=0.4
    a = 1/128.0 * (93.0 + 70 * alpha)
    b = 1.0/16.0 * (7.0 + 18.0 * alpha)
    c = 1.0/32.0 * (-7.0 + 14.0 *alpha)
    d = 1.0/16.0 * (1.0 - 2.0 * alpha)
    e = 1.0/128.0 * (-1.0 + 2.0 *alpha)

    Erow.fill(0.0); B8up.fill(0.0)

    Erow[-4] = e/2; Erow[-3] = d/2; Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a
    Erow[4] = e/2; Erow[3] = d/2; Erow[2] = c/2; Erow[1] = b/2

    a82 = np.array([1./2+alpha, (1./2+alpha)*2.0, 1./2+alpha, 0.0, 0.0, 0.0, 0.0])/2.
    a83 = np.array([-1/8.+alpha/4., 1./2+alpha, 2.0*(5./8+3*alpha/4.), 1./2+alpha, -1/8.+alpha/4., 0.0, 0.0])/2.
    a84 = np.array([1./32-alpha/16., -3./16+3*alpha/8., 15./32+17*alpha/16., 2.0*(11/16.+5.*alpha/8.), \
            15./32+17*alpha/16., -3./16+3*alpha/8., 1./32-alpha/16.])/2.

    for i in range(len(a84)):
        B8up[0, i] = a82[i]
        B8up[1, i] = a83[i]
        B8up[2, i] = a84[i]

    B8down = np.fliplr(np.flipud(B8up))
    explicit8loc = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit8 = np.zeros((pts, pts))

    for i in range(pts):
        explicit8loc[i, :] = np.roll(Erow, i)
        implicit8[i, :] = np.roll(Irow, i)

    explicit8loc[0, :]  = 0.0; explicit8loc[0, 0]   = 1.0
    explicit8loc[-1, :] = 0.0; explicit8loc[-1, -1] = 1.0
    implicit8[0, :]  = 0.0; implicit8[0, 0]   = 1.0
    implicit8[-1, :] = 0.0; implicit8[-1, -1] = 1.0
    explicit8loc[1:4, :] = B8up
    explicit8loc[-4:-1, :] = B8down

    A8loc = la.solve(implicit8, explicit8loc)
    eig8loc = la.eigvals(A8loc)

    eig8loc = np.sort(np.real(eig8loc))[::-1]
    eig8loc = (eig8loc - np.min(eig8loc))/np.ptp(eig8loc)
    imagloc = np.imag(eig8loc)

    print('8th order compact loc, max=%.15f, min=%.15f'%(np.max(np.real(eig8loc)), np.min(np.real(eig8loc))))

    # 4th order compact (DRP)
    a = 0.243527493120
    b = -0.204788880640
    c = 0.120007591680
    d = -0.045211119360
    e = 0.008228661760

    Erow.fill(0.0)

    Erow[-4] = e; Erow[-3] = d; Erow[-2] = c; Erow[-1] = b; Erow[0] = a
    Erow[4] = e; Erow[3] = d; Erow[2] = c; Erow[1] = b

    a42 = np.array([-0.057717512738,0.199278374994,-0.292668277650,0.244537361546,\
            -0.134605018019,0.056184263460,-0.015009191593, 0.0, 0.0, 0.0, 0.0])
    a43 = np.array([0.032649010764,-0.143339502575,0.273321177980,-0.294622121167,\
            0.186711738069,-0.062038376258,0.007318073189, 0.0, 0.0, 0.0, 0.0])
    a44 = np.array([-0.000054596010, 0.042124772446, -0.173103107841, 0.299615871352, -0.276543612935,\
    0.131223506571,-0.023424966418,0.013937561779,-0.024565095706,0.013098287852,-0.002308621090])

    B4up = np.zeros((3, pts))
    for i in range(len(a44)):
        B4up[0, i] = a42[i]
        B4up[1, i] = a43[i]
        B4up[2, i] = a44[i]

    B4down = np.fliplr(np.flipud(B4up))
    explicit4 = np.zeros((pts, pts))

    for i in range(pts):
        explicit4[i, :] = np.roll(Erow, i)

    explicit4[0, :]  = 0.0; explicit4[0, 0]   = 1.0
    explicit4[-1, :] = 0.0; explicit4[-1, -1] = 1.0
    explicit4[1:4, :] = B4up
    explicit4[-4:-1, :] = B4down

    eig4 = la.eigvals(np.eye(pts)-explicit4)

    imag = np.imag(eig4)[::-1]
    pos_imag_eig4 = imag[imag>=0]
    kexact_pos4 =kexact[np.where(imag>=0)]
    neg_imag_eig4 = imag[imag<0]
    kexact_neg4 =kexact[np.where(imag<0)]

    print('4th order optimized finite-difference, max=%.15f, min=%.15f'%(np.max(np.real(eig4)), np.min(np.real(eig4))))

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    # plot real
    plt.plot(kexact, np.sort(np.real(eig4))[::-1], color='mediumorchid', linestyle='--', linewidth=3, \
            label='4th order optimized finite-difference')
    plt.plot(kexact, (np.real(eig8loc)), 'g--',  linewidth=3, label='8th order compact (low order centered, $\\alpha=0.49$)')
    plt.plot(kexact, (np.real(eig8))[::-1], 'r--',  linewidth=3, label='8th order compact (8th order boundary, $\\alpha=0.49$)')
    plt.plot(kexactB, (np.real(eig)), 'b--',  linewidth=3,  label='8th order B-spline ($B\'\'\'\', \lambda^2 = 10^{-6}$)')

    plt.plot(kexact, np.ones(len(kexact)), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=np.min(np.real(eig8))-0.1, top=1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (real)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('tf-real-compact-8-bspline.eps', transparent=True, dpi=300)

    # plot imaginary
    plt.rc('text', usetex=True)
    plt.figure(2, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    plt.plot(kexact_pos4, pos_imag_eig4, color='mediumorchid', linestyle='--',\
            linewidth=3, alpha=1.0, label='4th order optimized finite-difference')
    plt.plot(kexact_neg4, neg_imag_eig4, '--',  linewidth=3, alpha=1.0)
    plt.plot(kexact, imagloc, 'g-', linewidth=4, alpha=1.0, label='8th order compact (low order centered, $\\alpha=0.49$)')
    plt.plot(kexact_pos, pos_imag_eig8, 'r--',  linewidth=3, alpha=1.0, label='8th order compact (8th order boundary, $\\alpha=0.49$)')
    plt.plot(kexact_neg, neg_imag_eig8, 'r--',  linewidth=3, alpha=1.0)
    plt.plot(kexactB, (np.imag(eig))[::-1], 'b:', linewidth=4, alpha=1.0, label='8th order B-spline ($B\'\'\'\', \lambda^2 = 10^{-6}$)')

    plt.ylim(bottom=-0.2, top=0.2)
    plt.xticks(fontsize=30)
    plt.yticks([-0.2, -0.1, 0.0, 0.10, 0.2],fontsize=30)
    plt.ylabel('Transfer function (imaginary)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.legend(fontsize=28, framealpha=0.0, loc='lower left')
    plt.tight_layout()
    plt.savefig('tf-imag-compact-8-bspline.pdf', transparent=True, dpi=300)

    # plot real-imaginary compact
    plt.rc('text', usetex=True)
    fig = plt.figure(3, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"
    ax = fig.add_subplot(111)

    plt.plot(np.real(eig4), np.imag(eig4), color='mediumorchid', marker='o', mfc='none', mew=2.0, \
            linestyle='none', label='4th order optimized finite-difference')
    ax.text(0.5, 0.2, 'Real: max = %0.6f, min = %0.6f'%(np.max(np.real(eig4)), np.min(np.real(eig4))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)
    ax.text(0.5, 0.1, 'Imag: max = %0.6f, min = %0.6f'%(np.max(np.imag(eig4)), np.min(np.imag(eig4))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)

    plt.axvline(x=0,  color='k', linestyle='--', linewidth=3)
    plt.axvline(x=1.0,  color='k', linestyle='--', linewidth=3)
    plt.ylim(bottom=-0.2, top=0.2)
    plt.xticks(fontsize=30)
    plt.yticks([-0.2, -0.1, 0.0, 0.10, 0.2],fontsize=30)
    plt.ylabel('Imaginary', labelpad=8, fontsize=30)
    plt.xlabel('Real', labelpad=8, fontsize=30)
    #plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('tf-ImRe-optimized-4.eps', transparent=True, dpi=300)

    # plot real-imaginary compact
    plt.rc('text', usetex=True)
    fig = plt.figure(4, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"
    ax = fig.add_subplot(111)

    plt.plot(np.real(eig8loc), np.imag(eig8loc), 'go', mfc='none', mew=2.0, label='8th order compact (low order centered, $\\alpha=0.49$)')
    ax.text(0.5, 0.2, 'Real: max = %0.6f, min = %0.6f'%(np.max(np.real(eig8loc)), np.min(np.real(eig8loc))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)
    ax.text(0.5, 0.1, 'Imag: max = %0.6f, min = %0.6f'%(np.max(np.imag(eig8loc)), np.min(np.imag(eig8loc))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)

    plt.axvline(x=0,  color='k', linestyle='--', linewidth=3)
    plt.axvline(x=1.0,  color='k', linestyle='--', linewidth=3)
    plt.ylim(bottom=-0.2, top=0.2)
    plt.xticks(fontsize=30)
    plt.yticks([-0.2, -0.1, 0.0, 0.10, 0.2],fontsize=30)
    plt.ylabel('Imaginary', labelpad=8, fontsize=30)
    plt.xlabel('Real', labelpad=8, fontsize=30)
    #plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('tf-ImRe-compact-8-loc.eps', transparent=True, dpi=300)

    # plot real-imaginary compact
    plt.rc('text', usetex=True)
    fig = plt.figure(5, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"
    ax = fig.add_subplot(111)

    plt.plot(np.real(eig8), np.imag(eig8), 'ro', mfc='none', mew=2.0, label='8th order compact (8th order boundary, $\\alpha=0.49$)')
    ax.text(0.5, 0.2, 'Real: max = %0.6f, min = %0.6f'%(np.max(np.real(eig8)), np.min(np.real(eig8))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)
    ax.text(0.5, 0.1, 'Imag: max = %0.6f, min = %0.6f'%(np.max(np.imag(eig8)), np.min(np.imag(eig8))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)

    plt.axvline(x=0,  color='k', linestyle='--', linewidth=3)
    plt.axvline(x=1.0,  color='k', linestyle='--', linewidth=3)
    plt.ylim(bottom=-0.2, top=0.2)
    plt.xticks(fontsize=30)
    plt.yticks([-0.2, -0.1, 0.0, 0.10, 0.2],fontsize=30)
    plt.ylabel('Imaginary', labelpad=8, fontsize=30)
    plt.xlabel('Real', labelpad=8, fontsize=30)
    #plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('tf-ImRe-compact-8.eps', transparent=True, dpi=300)

    # plot real-imaginary B-spline
    plt.rc('text', usetex=True)
    fig = plt.figure(6, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"
    ax = fig.add_subplot(111)

    plt.plot(np.real(eig), np.imag(eig), 'bo', mfc='none', mew=2.0, label='8th order B-spline ($B\'\'\'\', \lambda^2 = 10^{-6}$)')
    ax.text(0.5, 0.2, 'Real: max = %0.6f, min = %0.6f'%(np.max(np.real(eig)), np.min(np.real(eig))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)
    ax.text(0.5, 0.1, 'Imag: max = %0.6f, min = %0.6f'%(np.max(np.imag(eig)), np.min(np.imag(eig))), \
    fontsize=32, horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)

    plt.axvline(x=0,  color='k', linestyle='--', linewidth=3)
    plt.axvline(x=1.0,  color='k', linestyle='--', linewidth=3)
    plt.ylim(bottom=-0.2, top=0.2)
    plt.xticks(fontsize=30)
    plt.yticks([-0.2, -0.1, 0.0, 0.10, 0.2],fontsize=30)
    plt.ylabel('Imaginary', labelpad=8, fontsize=30)
    plt.xlabel('Real', labelpad=8, fontsize=30)
    #plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('tf-ImRe-bspline.eps', transparent=True, dpi=300)

    plt.show()
    return


def test_eig_periodic_transfer(pts):
    import numpy.linalg as la
    import scipy.linalg as sla

    alpha = 0.49
    dx = 2.0*np.pi/(pts-1)
    kexact = np.linspace(0.0, np.pi, pts)

    # 4th order compact
    a = 5/8.0 + 0.75 * alpha
    b = 1.0/2.0 +  alpha
    c = -1.0/8.0 + 0.25 * alpha

    Erow = np.zeros(pts)
    Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a; Erow[1] = b/2; Erow[2] = c/2
    explicit4 = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit4 = np.zeros((pts, pts))

    for i in range(pts):
        explicit4[i, :] = np.roll(Erow, i)
        implicit4[i, :] = np.roll(Irow, i)

    A4 = la.solve(implicit4, explicit4)
    eig4 = la.eigvals(A4)

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

    A8 = la.solve(implicit8, explicit8)
    eig8 = la.eigvals(A8)

    # plot real
    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    plt.plot(kexact, np.sort(np.real(eig4))[::-1], 'b--', linewidth=3, label='4th order compact ($\\alpha=0.49$)')
    plt.plot(kexact, np.sort(np.imag(eig4))[::-1], 'b:',  linewidth=3, label='4th order compact ($\\alpha=0.49$)')
    plt.plot(kexact, np.sort(np.real(eig8))[::-1], 'r--', linewidth=3, label='8th order compact ($\\alpha=0.49$)')
    plt.plot(kexact, np.sort(np.imag(eig8))[::-1], 'r:',  linewidth=3, label='8th order compact ($\\alpha=0.49$)')

    plt.ylim(bottom=-0.1, top=1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    #plt.savefig('transfer-compact-8-area-bspline.eps', transparent=True, dpi=300)
    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')

    # plot real-imaginary
    plt.rc('text', usetex=True)
    plt.figure(2, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    plt.plot(np.real(eig4), np.imag(eig4), 'bo', mfc='none', mew=1.0, linewidth=3, label='4th order compact ($\\alpha=0.49$)')
    plt.plot(np.real(eig8), np.imag(eig8), 'ro', mfc='none', mew=1.0, linewidth=3, label='8th order compact ($\\alpha=0.49$)')

    plt.ylim(bottom=-0.1, top=0.1)
    plt.xticks(fontsize=30)
    plt.yticks([-0.1, 0.0, 0.1],fontsize=30)
    plt.ylabel('Imaginary', labelpad=8, fontsize=30)
    plt.xlabel('Real', labelpad=8, fontsize=30)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    #plt.savefig('transfer-compact-8-area-bspline.eps', transparent=True, dpi=300)

    plt.show()

    return


def test_smooth_spline_vs_optimized_fd_transfer_area(degree, pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    degree=10
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
#        _ ,_ ,_ ,_ , fullP_1D = \
#        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'periodic', matrix_type='full')

        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)
    else:
        _ ,_ ,_ ,_ , fullP_1D = \
        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)

    lam = np.array([1e-5])

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)
    yplot = np.zeros((3, median))
    yplot_compact = np.zeros((2, pts-1))

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

#   optimized finite-difference filter coefficients

    compactlabel=['Optimized (9pt stencil)','4th order optimized (13 point stencil)']
    color=['r','r--']
    for j in range(1):
        if j==0:
            d0 = 0.243527493120
            d1 = -0.204788880640
            d2 = 0.120007591680
            d3 = -0.045211119360
            d4 = 0.008228661760
            d5 = 0.0
            d6 = 0.0
        elif j==1:
            d0 = 0.190899511506
            d1 = -0.171503832236
            d2 = 0.123632891797
            d3 = -0.069975429105
            d4 = 0.029662754736
            d5 = -0.008520738659
            d6 = 0.001254597714

        for i in range(pts):
            transfer[i] = 1.0 - (d0 + 2 * d1 * np.cos(kexact[i]) + 2*d2 * np.cos(2.0 * kexact[i]) + \
                2*d3 * np.cos(3.0 * kexact[i]) + 2*d4 * np.cos(4.0 * kexact[i]) + \
                2*d5 * np.cos(5.0 * kexact[i]) + 2*d6 * np.cos(6.0 * kexact[i]))

#        plt.plot(kexact[:-1], transfer[:-1], 'k--', linewidth=3, label=compactlabel[j])
        plt.plot(kexact, transfer, color[j], linewidth=3, label=compactlabel[j])
        yplot_compact[j, :] = transfer[:-1]

    lam = np.array([1e-1, 1e-13, 2e-3])
    lamlabel=['10^{-1}', '10^{-6}', '5 \\times10^{-3}']
    color=['b:','b--', 'b-.']
    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], color[j], linewidth=3, label='8th-order B-spline \
                ($B^{(6)}, \lambda^2 = %s)$'%lamlabel[j])
        yplot[j, :] = transfer[:median]

    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('$k\Delta x$', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.fill_between(x[:median], yplot[0,:], yplot[1,:], color='b', alpha=0.2)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('transfer-optimized-4-area-bspline.pdf', transparent=True, dpi=300)
    plt.show()
    return


def test_smooth_spline_vs_fd_transfer_area(degree, pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
#        _ ,_ ,_ ,_ , fullP_1D = \
#        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'periodic', matrix_type='full')

        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)
    else:
        _ ,_ ,_ ,_ , fullP_1D = \
        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)
    yplot = np.zeros((2, median))
    yplot_compact = np.zeros((2, pts-1))

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

#    compact filter coefficients

    a1 = np.array([0.0, 0.49])
    b1 = np.array([0.0, 0.0])
    compactlabel=['4th-order explicit ($\\alpha=0$)',\
                  '4th-order compact ($\\alpha=0.49$)']
    color=['r:', 'r--']
    for j in range(len(a1)):
        if j==0 or j==2:
            a = 5/8.0 + 0.75 * a1[j]
            b = 1.0/2.0 +  a1[j]
            c = -1.0/8.0 + 0.25 * a1[j]
            d = 0.0
            e = 0.0
        elif j==1 or j==3:
            a = 1/128.0 * (93.0 + 70 * a1[j])
            b = 1.0/16.0 * (7.0 + 18.0 * a1[j])
            c = 1.0/32.0 * (-7.0 + 14.0 * a1[j])
            d = 1.0/16.0 * (1.0 - 2.0 * a1[j])
            e = 1.0/128.0 * (-1.0 + 2.0 * a1[j])

        num = 0.0; den = 0.0
        for i in range(pts):
            num = a + b * np.cos(kexact[i]) + c * np.cos(2.0 * kexact[i]) + \
                d * np.cos(3.0 * kexact[i]) + e * np.cos(4.0 * kexact[i])

            den = 1.0 + 2.0 * a1[j] * np.cos(kexact[i]) + 2.0 * b1[j] * np.cos(2.0 * kexact[i])

            transfer[i] = num / den

#        plt.plot(kexact[:-1], transfer[:-1], 'k--', linewidth=3, label=compactlabel[j])
        plt.plot(kexact[:-1], transfer[:-1], color[j], linewidth=3, label=compactlabel[j])
        yplot_compact[j, :] = transfer[:-1]

    plt.fill_between(kexact[:-1], yplot_compact[0,:], yplot_compact[1,:], color='r', alpha=0.2)

    lam = np.array([1e-1, 1e-13])
    lamlabel=['10^{-1}', '10^{-6}']
    color=['b:', 'b--']
    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], color[j], linewidth=3, label='4th-order B-spline ($B\'\', \lambda^2 = %s)$'%lamlabel[j])
        yplot[j, :] = transfer[:median]

    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('$k\Delta x$', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.fill_between(x[:median], yplot[0,:], yplot[1,:], color='b', alpha=0.2)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
#    plt.savefig('transfer-compact-4-area-bspline.pdf', transparent=True, dpi=300)


    plt.rc('text', usetex=True)
    plt.figure(2, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

#    compact filter coefficients

    yplot = np.zeros((2, median))
    yplot_compact = np.zeros((2, pts-1))

    a1 = np.array([0.0, 0.49])
    b1 = np.array([0.0, 0.0])
    compactlabel=['8th-order explicit ($\\alpha=0$)',\
                  '8th-order compact ($\\alpha=0.49$)']
    color=['r:', 'r--']
    for j in range(len(a1)):
        if j==0 or j==2:
            a = 5/8.0 + 0.75 * a1[j]
            b = 1.0/2.0 +  a1[j]
            c = -1.0/8.0 + 0.25 * a1[j]
            d = 0.0
            e = 0.0
        elif j==1 or j==3:
            a = 1/128.0 * (93.0 + 70 * a1[j])
            b = 1.0/16.0 * (7.0 + 18.0 * a1[j])
            c = 1.0/32.0 * (-7.0 + 14.0 * a1[j])
            d = 1.0/16.0 * (1.0 - 2.0 * a1[j])
            e = 1.0/128.0 * (-1.0 + 2.0 * a1[j])

        num = 0.0; den = 0.0
        for i in range(pts):
            num = a + b * np.cos(kexact[i]) + c * np.cos(2.0 * kexact[i]) + \
                d * np.cos(3.0 * kexact[i]) + e * np.cos(4.0 * kexact[i])

            den = 1.0 + 2.0 * a1[j] * np.cos(kexact[i]) + 2.0 * b1[j] * np.cos(2.0 * kexact[i])

            transfer[i] = num / den

#        plt.plot(kexact[:-1], transfer[:-1], 'k--', linewidth=3, label=compactlabel[j])
        plt.plot(kexact[:-1], transfer[:-1], color[j], linewidth=3, label=compactlabel[j])
        yplot_compact[j, :] = transfer[:-1]

    plt.fill_between(kexact[:-1], yplot_compact[0,:], yplot_compact[1,:], color='r', alpha=0.2)

    degree = 10
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)

    lam = np.array([1e-1, 1e-13])
    lamlabel=['10^{-1}', '10^{-6}']
    color=['b:', 'b--']
    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], color[j], linewidth=3, label='8th-order B-spline \
                ($B^{(6)}, \lambda^2 = %s)$'%lamlabel[j])
        yplot[j,:] = transfer[:median]


    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('$k\Delta x$', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.fill_between(x[:median], yplot[0,:], yplot[1,:], color='b', alpha=0.2)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('transfer-compact-8-area-bspline.pdf', transparent=True, dpi=300)
    plt.show()
    return


def test_smooth_spline_vs_fd_transfer(degree, pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
#        _ ,_ ,_ ,_ , fullP_1D = \
#        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'periodic', matrix_type='full')

        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)
    else:
        _ ,_ ,_ ,_ , fullP_1D = \
        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)

    lam = np.array([1e-5])

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

#    compact filter coefficients

    a1 = np.array([0.0, 0.0, 0.49, 0.49])
    b1 = np.array([0.0, 0.0, 0.0, 0.0])
    compactlabel=['4th order explicit ($\\alpha=0$)',\
                  '8th order explicit ($\\alpha=0$)',\
                  '4th order compact ($\\alpha=0.49$)',\
                  '8th order compact ($\\alpha=0.49$)']
    color=['b:', 'r:', 'b--', 'r--']
#    a1 = np.array([0.0, 0.49])
#    b1 = np.array([0.0, 0.0])
#    compactlabel=['8th order explicit ($\\alpha=0$)',\
#                  '8th order compact ($\\alpha=0.49$)']
#    color=['r:', 'r--']
    for j in range(len(a1)):
        if j==0 or j==2:
            a = 5/8.0 + 0.75 * a1[j]
            b = 1.0/2.0 +  a1[j]
            c = -1.0/8.0 + 0.25 * a1[j]
            d = 0.0
            e = 0.0
        elif j==1 or j==3:
            a = 1/128.0 * (93.0 + 70 * a1[j])
            b = 1.0/16.0 * (7.0 + 18.0 * a1[j])
            c = 1.0/32.0 * (-7.0 + 14.0 * a1[j])
            d = 1.0/16.0 * (1.0 - 2.0 * a1[j])
            e = 1.0/128.0 * (-1.0 + 2.0 * a1[j])

        num = 0.0; den = 0.0
        for i in range(pts):
            num = a + b * np.cos(kexact[i]) + c * np.cos(2.0 * kexact[i]) + \
                d * np.cos(3.0 * kexact[i]) + e * np.cos(4.0 * kexact[i])

            den = 1.0 + 2.0 * a1[j] * np.cos(kexact[i]) + 2.0 * b1[j] * np.cos(2.0 * kexact[i])

            transfer[i] = num / den

#        plt.plot(kexact[:-1], transfer[:-1], 'k--', linewidth=3, label=compactlabel[j])
        plt.plot(kexact[:-1], transfer[:-1], color[j], linewidth=3, label=compactlabel[j])

    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], 'b', linewidth=3, label='4th order B-spline ($B\'\', \lambda^2 = 10^{-4})$')

    degree = 10
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)

    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], 'r', linewidth=3, label='8th order B-spline ($B\'\'\'\', \lambda^2 = 10^{-4})$')

#    plt.annotate(r'$a$', (1.2, 0.7), fontsize=35)
#    plt.annotate(r'$b$', (1.85, 0.7), fontsize=35)
#    plt.annotate(r'$c$', (2.3, 0.7), fontsize=35)
#    plt.annotate(r'$d$', (2.52, 0.7), fontsize=35)

    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.tight_layout()
    plt.legend(fontsize=28, framealpha=0.0)
    plt.savefig('transfer-compact-bspline.eps', transparent=True, dpi=300)
    plt.show()
    return

def test_smooth_spline_vs_drp7_optimized_transfer(degree, pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
#        _ ,_ ,_ ,_ , fullP_1D = \
#        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'periodic', matrix_type='full')

        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)
    else:
        _ ,_ ,_ ,_ , fullP_1D = \
        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)

    lam = np.array([1e-5])

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

#   seven-point DRP filter coefficients

    #compactlabel=['4th order (compact)','8th order (compact)']
    compactlabel=['DRP seven-point stencil 1','DRP seven-point stencil 2']
    color=['b:','b--']
    for j in range(2):
        if j==0:
            d0 = 0.327698660846
            d1 = -0.235718815308
            d2 = 0.086150669577
            d3 = -0.014281184692
        elif j==1:
            d0 = 0.303058855153
            d1 = -0.231695528614
            d2 = 0.098470572424
            d3 = -0.018304471386

        for i in range(pts):
            transfer[i] = 1.0 - (d0 + 2 * d1 * np.cos(kexact[i]) + 2*d2 * np.cos(2.0 * kexact[i]) + \
                2*d3 * np.cos(3.0 * kexact[i]))

#        plt.plot(kexact[:-1], transfer[:-1], 'k--', linewidth=3, label=compactlabel[j])
        plt.plot(kexact, transfer, color[j], linewidth=3, label=compactlabel[j])

    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], 'b', linewidth=3, label='4th order B-spline ($B\'\', \lambda^2 = 10^{-4}$)')

#    degree = 9
#    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)
#
#    fullB_1D = ders[:, 0, :]
#    fullI_1D = np.eye(fullB_1D.shape[0])
#
#    if periodic == True:
#        row, col = fullB_1D.shape
#        periodic_row = np.zeros(row)
#        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]
#
#        periodic_row = np.roll(periodic_row, -2)
#
#        for i in range(row):
#            fullB_1D[i, :] = np.roll(periodic_row, i)
#
#    for j in range(len(lam)):
#        for i in range(pts):
#            input = np.exp(i* 1j* x)
#            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))
#
#            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
#            transfer[i] = np.mean(np.abs(output)/np.abs(input))
#
#        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
##        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
#        plt.plot(x[:median], transfer[:median], 'r', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.0e'%(degree-1, lam[j]*10))
#
#
##    plt.annotate(r'$a$', (1.2, 0.7), fontsize=35)
##    plt.annotate(r'$b$', (1.85, 0.7), fontsize=35)
##    plt.annotate(r'$c$', (2.3, 0.7), fontsize=35)
##    plt.annotate(r'$d$', (2.52, 0.7), fontsize=35)

    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.tight_layout()
    plt.legend(fontsize=28, framealpha=0.0)
    plt.savefig('transfer-optimized-drp7-bspline.eps', transparent=True, dpi=300)
    plt.show()
    return



def test_smooth_spline_vs_optimized_fd_transfer(degree, pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
#        _ ,_ ,_ ,_ , fullP_1D = \
#        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'periodic', matrix_type='full')

        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)
    else:
        _ ,_ ,_ ,_ , fullP_1D = \
        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)

    lam = np.array([1e-5])

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

#   optimized finite-difference filter coefficients

    #compactlabel=['4th order (compact)','8th order (compact)']
    compactlabel=['4th order optimized (9 point stencil)','4th order optimized (13 point stencil)']
    color=['b:','b--']
    for j in range(2):
        if j==0:
            d0 = 0.243527493120
            d1 = -0.204788880640
            d2 = 0.120007591680
            d3 = -0.045211119360
            d4 = 0.008228661760
            d5 = 0.0
            d6 = 0.0
        elif j==1:
            d0 = 0.190899511506
            d1 = -0.171503832236
            d2 = 0.123632891797
            d3 = -0.069975429105
            d4 = 0.029662754736
            d5 = -0.008520738659
            d6 = 0.001254597714

        for i in range(pts):
            transfer[i] = 1.0 - (d0 + 2 * d1 * np.cos(kexact[i]) + 2*d2 * np.cos(2.0 * kexact[i]) + \
                2*d3 * np.cos(3.0 * kexact[i]) + 2*d4 * np.cos(4.0 * kexact[i]) + \
                2*d5 * np.cos(5.0 * kexact[i]) + 2*d6 * np.cos(6.0 * kexact[i]))

#        plt.plot(kexact[:-1], transfer[:-1], 'k--', linewidth=3, label=compactlabel[j])
        plt.plot(kexact, transfer, color[j], linewidth=3, label=compactlabel[j])

    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], 'b', linewidth=3, label='4th order B-spline ($B\'\', \lambda^2 = 10^{-4}$)')

#    degree = 9
#    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)
#
#    fullB_1D = ders[:, 0, :]
#    fullI_1D = np.eye(fullB_1D.shape[0])
#
#    if periodic == True:
#        row, col = fullB_1D.shape
#        periodic_row = np.zeros(row)
#        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]
#
#        periodic_row = np.roll(periodic_row, -2)
#
#        for i in range(row):
#            fullB_1D[i, :] = np.roll(periodic_row, i)
#
#    for j in range(len(lam)):
#        for i in range(pts):
#            input = np.exp(i* 1j* x)
#            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))
#
#            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
#            transfer[i] = np.mean(np.abs(output)/np.abs(input))
#
#        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
##        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
#        plt.plot(x[:median], transfer[:median], 'r', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.0e'%(degree-1, lam[j]*10))
#
#
##    plt.annotate(r'$a$', (1.2, 0.7), fontsize=35)
##    plt.annotate(r'$b$', (1.85, 0.7), fontsize=35)
##    plt.annotate(r'$c$', (2.3, 0.7), fontsize=35)
##    plt.annotate(r'$d$', (2.52, 0.7), fontsize=35)

    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.tight_layout()
    plt.legend(fontsize=28, framealpha=0.0)
    plt.savefig('transfer-optimized-bspline.eps', transparent=True, dpi=300)
    plt.show()
    return


def test_smooth_spline_vs_compact_transfer(degree, pts, periodic=True):
    import pygsvd
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
#        _ ,_ ,_ ,_ , fullP_1D = \
#        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'periodic', matrix_type='full')

        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)
    else:
        _ ,_ ,_ ,_ , fullP_1D = \
        smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)

    lam = np.array([1e-5])

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

#    compact filter coefficients

    a1 = np.array([0.35, 0.45])
    b1 = np.array([0.0, 0.0])
    #compactlabel=['4th order (compact)','8th order (compact)']
    compactlabel=['4th order (Finite-difference)','8th order (Finite-difference)']
    color=['b--','r--']
    for j in range(2):
        if j==0:
            a = 5/8.0 + 0.75 * a1[j]
            b = 1.0/2.0 +  a1[j]
            c = -1.0/8.0 + 0.25 * a1[j]
            d = 0.0
            e = 0.0
        elif j==1:
            a = 1/128.0 * (93.0 + 70 * a1[j])
            b = 1.0/16.0 * (7.0 + 18.0 * a1[j])
            c = 1.0/32.0 * (-7.0 + 14.0 * a1[j])
            d = 1.0/16.0 * (1.0 - 2.0 * a1[j])
            e = 1.0/128.0 * (-1.0 + 2.0 * a1[j])

        num = 0.0; den = 0.0
        for i in range(pts):
            num = a + b * np.cos(kexact[i]) + c * np.cos(2.0 * kexact[i]) + \
                d * np.cos(3.0 * kexact[i]) + e * np.cos(4.0 * kexact[i])

            den = 1.0 + 2.0 * a1[j] * np.cos(kexact[i]) + 2.0 * b1[j] * np.cos(2.0 * kexact[i])

            transfer[i] = num / den

#        plt.plot(kexact[:-1], transfer[:-1], 'k--', linewidth=3, label=compactlabel[j])
        plt.plot(kexact[:-1], transfer[:-1], color[j], linewidth=3, label=compactlabel[j])

    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], 'b', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.0e'%(degree-1, lam[j]*10))

    degree = 9
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)

    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))
        plt.plot(x[:median], transfer[:median], 'r', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.0e'%(degree-1, lam[j]*10))


#    plt.annotate(r'$a$', (1.2, 0.7), fontsize=35)
#    plt.annotate(r'$b$', (1.85, 0.7), fontsize=35)
#    plt.annotate(r'$c$', (2.3, 0.7), fontsize=35)
#    plt.annotate(r'$d$', (2.52, 0.7), fontsize=35)

    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Exact')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.tight_layout()
    plt.legend(fontsize=28)
    #plt.savefig('transfer-compact-bspline.eps', transparent=True, dpi=300)
    plt.show()
    return

def test_smooth_spline_transfer_function_derivative(pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)


    lam = np.array([1e-3, 2.5e-4, 1e-16])
    degreelist = np.array([8, 9, 11])
    #lam = np.array([1e-3, 1e-8, 1e-16])
    #degreelist = np.array([6, 9, 15])


    derivlabel = np.array(['B^{(2)}',\
                           'B^{(4)}',\
                           'B^{(6)}'])

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    color=['b','g','r']
    for k in range(len(degreelist)):
        degree = int(degreelist[k])
        fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

        fullB_1D = ders[:, 0, :]
        fullI_1D = np.eye(fullB_1D.shape[0])

        if periodic == True:
            row, col = fullB_1D.shape
            periodic_row = np.zeros(row)
            periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

            periodic_row = np.roll(periodic_row, -2)

            for i in range(row):
                fullB_1D[i, :] = np.roll(periodic_row, i)

        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[k] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
        plt.plot(x[:median], transfer[:median], color=color[k], linewidth=3, label='8th-order B-spline ($%s, \lambda^2 = 10^{-6}$)'%(derivlabel[k]))

    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('$k\Delta x$', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.tight_layout()
    plt.legend(fontsize=28, framealpha=0.0)
    plt.savefig('transfer-deriv.eps', transparent=True, dpi=300)
    plt.show()
    return


def test_smooth_spline_transfer_function_lambda(pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    x = np.linspace(mink, maxk, pts, dtype=np.float128)
    kexact = np.linspace(0.0, np.pi, pts)
    transfer = np.linspace(mink, maxk, pts,dtype=np.float128)

    lam = np.array([1e-1, 1e-2, 1e-3, 1e-4])
    lamlabel = np.array(['10^{-1}','10^{-2}','10^{-3}','10^{-4}'])

    dx = x[1] - x[0]
    median = int(pts/2)
    transfer = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8.1))
    plt.rcParams["font.family"] = "Times New Roman"

    degree = 10
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    if periodic == True:
        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)

    for j in range(len(lam)):
        for i in range(pts):
            input = np.exp(i* 1j* x)
            output = fullB_1D @ (sla.solve(fullB_1D.T @ fullB_1D + lam[j] * fullI_1D, fullB_1D.T @ input))

            #transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
            transfer[i] = np.mean(np.abs(output)/np.abs(input))

        transfer = (transfer - np.min(transfer)) / np.ptp(transfer)
#        plt.plot(x[:median], transfer[:median], 'k', linewidth=3, label='%dth order (Bspline), $\lambda^2$ = %.1e'%(degree-1, lam[j]))

        plt.plot(x[:median], transfer[:median], linewidth=3, label='8th-order B-spline \
                ($B^{(6)}, \lambda^2 = %s$)'%(lamlabel[j]))


#    plt.annotate(r'$a$', (1.2, 0.7), fontsize=35)
#    plt.annotate(r'$b$', (1.85, 0.7), fontsize=35)
#    plt.annotate(r'$c$', (2.3, 0.7), fontsize=35)
#    plt.annotate(r'$d$', (2.52, 0.7), fontsize=35)
    plt.plot(kexact, np.ones(pts), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('$k\Delta x$', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.tight_layout()
    plt.legend(fontsize=27, framealpha=0.0)
    plt.savefig('transfer-lambda.eps', transparent=True, dpi=300)
    plt.show()
    return

def test_smooth_spline_svd_coeff(degree, pts, periodic=False):
    import pygsvd
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    _ ,_ ,_ ,_ , fullP_1D = \
    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    if periodic == True:
        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)

    U, S, VT = la.svd(fullB_1D)
    C2, S2, X, U2, V2= pygsvd.gsvd(fullB_1D, fullP_1D, extras='uv', full_matrices=True)

    b = np.sin(4.0 * fullu)
    x = la.solve(fullB_1D, b)

    V = VT.T
    num = np.sort(abs(VT @ x))[::-1]
    num2 = np.sort(abs(la.solve(X, x)))[::-1]

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(12, 8))
    plt.semilogy(np.linspace(0, 19, 20, dtype=int), num[:20], \
            'o-', linewidth=3, mfc='none', ms=10, mew=2, label='SVD basis $|\chi_i|$')
    plt.semilogy(np.linspace(0, 19, 20, dtype=int), num2[:20],\
            'o--', linewidth=3, mfc='none', ms=10, mew=2, label='GSVD basis $|\chi_i\'|$')
    plt.xticks(np.arange(0, 20, step=5), fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(fontsize=40)

    a = plt.axes([.65, .3, .2, .2], facecolor='y')
    plt.plot(fullu, b, linewidth=3)
    plt.title('Exact Solution', fontsize=25)
    plt.xticks(np.arange(mink, maxk, step=2.0), fontsize=20)
    plt.yticks(np.arange(-1,1.01, step=1), fontsize=20)

    plt.tight_layout()
    plt.savefig('gsvd-svd-coeff.eps', transparent=True, dpi=300)
    plt.show()
    return

def test_smooth_spline_right_vector(degree, pts, periodic=False):
    import pygsvd
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    _ ,_ ,_ ,_ , fullP_1D = \
    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    if periodic == True:
        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)

    U, S, VT = la.svd(fullB_1D)
    C2, S2, X, U2, V2= pygsvd.gsvd(fullB_1D, fullP_1D, extras='uv', full_matrices=True)

    V = VT.T

    b = np.zeros(len(fullu))

    for i in range(len(fullu)):
        if fullu[i] > np.pi/2 and fullu[i] < np.pi:
            b[i] = np.sin(4.0 * fullu[i]) + np.random.random(1)
        else:
            b[i] = np.sin(4.0 * fullu[i])

    K = np.linspace(25, pts, 4, dtype=int)
    x = np.zeros((len(K),len(fullu)))
    x2 = np.zeros((len(K),len(fullu)))

    plt.rc('text', usetex=True)
    for i in range(len(K)):
        for k in range(K[i]):
            x[i, :] = x[i, :] + np.dot(U[:, k], b) / S[k] * V[:, k]
            x2[i, :] = x2[i, :] + np.dot(U2[:, k], b) / C2[k] * X[:, k] * (C2[k])**2

        plt.figure(i+1, figsize=(16, 14))
        plt.plot(fullu, x[i, :], '--', color='magenta', linewidth=3.5, label="SVD basis")
        plt.plot(fullu, x2[i, :], '--', color='mediumblue', linewidth=3.5, label="GSVD basis")
        plt.plot(fullu, b, 'black', linewidth=2.5, label="Noisy data")
        plt.title("$K = %d$"%K[i], fontsize=50)
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.ylim(top=1.8, bottom=-1.8)
        plt.legend(loc="lower right", fontsize=40)
        plt.tight_layout()
        plt.savefig('gsvd-svd-%d.eps'%K[i], transparent=True, dpi=300)
    plt.show()
    return

def test_smooth_spline_svd(degree, pts, periodic=False, penalty=False):
    import pygsvd
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 1.0
    # number of derivatives
    deriv = 2

    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)

    fullB_1D = ders[:, 0, :]
    fullI_1D = np.eye(fullB_1D.shape[0])

    _ ,_ ,_ ,_ , fullP_1D = \
    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, fullB_1D, 2, 'nonperiodic', matrix_type='full')

    if periodic == True:
        row, col = fullB_1D.shape
        periodic_row = np.zeros(row)
        periodic_row[:degree] = fullB_1D[degree+1, degree-1:2*degree-1]

        periodic_row = np.roll(periodic_row, -2)

        for i in range(row):
            fullB_1D[i, :] = np.roll(periodic_row, i)

    if penalty==False:
        U, S, VT = la.svd(fullB_1D)

        idx = np.argsort(abs(S))[::-1]
        U = U[:, idx]
        S = S[idx]
        VT = VT[:,idx]

        V = VT.T

    elif penalty==True:
        U, S, V, U2, V2= pygsvd.gsvd(fullB_1D, fullP_1D, extras='uv', full_matrices=True)
        VT = V.T
        V = VT.T * U
        S = U.copy()

    x = np.linspace(0, pts, len(S))

    plt.rc('text', usetex=True)
    fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(20, 10))
    for i in range(10):
        plt.subplot(2, 5, i+1)
        plt.plot(x, np.real(V[:,i]))
        plt.ylim(top=0.2, bottom=-0.2)
        plt.xticks(np.arange(0, pts, step=50), fontsize=15)
        plt.yticks(np.arange(-0.2, 0.22, step=0.2), fontsize=15)

        if penalty==False:
            plt.title('$i = {%d}, ~\sigma_{%d} = %0.4f$'%(i+1, i+1, S[i].real), \
                    fontsize=20, bbox=dict(facecolor='none', edgecolor='black', pad=5.0), position=(0.58, 1.02))
        elif penalty==True:
            plt.title('$i = {%d}, ~\sigma_{%d} = %0.2e$'%(i+1, i+1, S[i].real), \
                    fontsize=20, bbox=dict(facecolor='none', edgecolor='black', pad=5.0), position=(0.58, 1.02))

    if penalty==False:
        plt.savefig('allv.eps', transparent=True, dpi=300)
    elif penalty==True:
        plt.savefig('allx.eps', transparent=True, dpi=300)

    fig1, ax1 = plt.subplots(nrows=2, ncols=5, figsize=(20, 10))
    for i in range(10):
        plt.subplot(2, 5, i+1)
        plt.plot(x, np.real(V[:,-(i+1)]), \
                label='$i = %d$, $\sigma = %.4f$'%(V.shape[0]-i, S[-(i+1)].real))
        plt.ylim(top=0.2, bottom=-0.2)
        plt.xticks(np.arange(0, pts, step=50), fontsize=15)
        plt.yticks(np.arange(-0.2, 0.22, step=0.2), fontsize=15)

        if penalty==False:
            plt.title('$i = {%d}, ~\sigma_{%d} = %0.4f$'%(V.shape[0]-i, V.shape[0]-i, S[-(i+1)].real), \
                    fontsize=20, bbox=dict(facecolor='none', edgecolor='black', pad=5.0), position=(0.56, 1.02))
        elif penalty==True:
            plt.title('$i = {%d}, ~\sigma_{%d} = %0.2e$'%(V.shape[0]-i, V.shape[0]-i, S[-(i+1)].real), \
                    fontsize=20, bbox=dict(facecolor='none', edgecolor='black', pad=5.0), position=(0.56, 1.02))

    if penalty==False:
        plt.savefig('allvbehind.eps', transparent=True, dpi=300)
    elif penalty==True:
        plt.savefig('allxbehind.eps', transparent=True, dpi=300)
    plt.show()
    return

def test_tikhonov_filter(degree, pts):

    lamda = np.linspace(-1, 1, pts)
    sigma = np.linspace(0, 30, pts)

    plt.rc('text', usetex=True)
    fig = plt.figure(1, figsize=(8, 6))
    plt.rc('text.latex', preamble=r"\usepackage[dvipsnames]{xcolor}")
    plt.rcParams["font.family"] = "Times New Roman"
    ax = fig.add_subplot(111)

    plt.loglog(sigma, sigma**2/(sigma**2 + lamda**2), linewidth=3)
    #plt.plot(0.95, 0.01, 'r|', ms=20, mew=2)
    #plt.plot(5, 0.01, 'r|', ms=20, mew=2)
    plt.axvline(x=0.95, color='r', linestyle='--', linewidth=2)
    plt.axvline(x=5, color='r', linestyle='--', linewidth=2)
    plt.annotate('\\textcolor{red}{$\sigma_i = \lambda$}', xy=(0.95, 0.01), xytext=(0.7, 0.1), fontsize=22, rotation=90)
    plt.annotate('\\textcolor{red}{$\sigma_i = 10 \lambda$}', xy=(5, 0.01), xytext=(3.8, 0.13), fontsize=22, rotation=90)
    plt.ylim(0.01, 1*1.2)
    plt.xlim(0.1, 30)
    plt.yticks(fontsize=20)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.tick_params(which='both', length=6, width=1)
    plt.xlabel('$\sigma_i$', fontsize=25)
    plt.ylabel('$\\varphi^{[\lambda]} = \sigma^2_i/ (\sigma^2_i + \lambda^2)$', fontsize=25)
    #plt.title('Filter Factors for Tikhonov regularization', fontsize=25)
    plt.tight_layout()
    plt.savefig('filter-factor.eps', transparent=True, dpi=300)
    plt.show()

    return

def test_approximation_curve():
    from curve import nurbs
    import numpy.linalg as la
    from findspan import findspan
    from basisfunctions import basisfunctions
    from knot_multiplicity import knot_multiplicity

    # weight
    w = 1.0
    # number of control pts
    n_cp = 8
    # degree
    degree = 5
    # number of pts
    pts = n_cp - degree + 1
    # knot vector
    mink = 0.0
    maxk = 1.0
    U = np.linspace(mink, maxk, pts)
    # knot multiplicity
    U = knot_multiplicity(degree, U)
    # dimensions
    dim = 3
    # n_splines is the number of B-spline curves
    n_splines = len(U) - degree - 1
    # points
    pts = n_splines
    # control points
    P = np.zeros((dim, n_cp))
    P[-1, :] = 1.0
    # interpolation points
    D = np.zeros((dim, n_cp))
    # curve points
    C = np.zeros((dim, pts))
    # basis functions
    N = np.zeros((n_splines, n_splines))

    assert(len(U) == (n_cp + degree + 1))

#    # initialize P_i
#    D[0:dim-1, 0] = np.array([0.1, 0.0])
#    D[0:dim-1, 1] = np.array([0.12, 0.25])
#    D[0:dim-1, 2] = np.array([0.02, 0.45])
#    D[0:dim-1, 3] = np.array([0.0, 0.8])
#    D[0:dim-1, 4] = np.array([0.3, 1.0])
#    D[0:dim-1, 5] = np.array([0.85, 0.95])
#    D[0:dim-1, 6] = np.array([1.1, 0.6])
#    D[0:dim-1, 7] = np.array([1.06, 0.3])
#    D[0:dim-1, 8] = np.array([0.95, 0.09])
#    D[0:dim-1, 9] = np.array([0.7, -0.08])

    D[0:dim-1, 0] = np.array([0.1, 0.0])
    D[0:dim-1, 1] = np.array([0.12, 0.25])
    D[0:dim-1, 2] = np.array([0.02, 0.45])
    D[0:dim-1, 3] = np.array([0.0, 0.8])
    D[0:dim-1, 4] = np.array([1.1, 0.6])
    D[0:dim-1, 5] = np.array([1.06, 0.3])
    D[0:dim-1, 6] = np.array([0.95, 0.09])
    D[0:dim-1, 7] = np.array([0.7, -0.08])

    bc0 = D[0:dim-1, 0]
    bc1 = D[0:dim-1, -1]

    totaldist = 0.0
    for j in range(1, n_cp):
        totaldist += np.sqrt((D[1, j-1] - D[1, j]) ** 2 + (D[0, j-1] - D[0, j]) ** 2)

    u = np.zeros(n_cp)
    u[0] = mink
    for j in range(1, len(u)):
        dist = np.sqrt((D[1, j-1] - D[1, j]) ** 2 + (D[0, j-1] - D[0, j]) ** 2) / totaldist
        u[j] = u[j-1] + dist

#    u = np.linspace(mink, maxk, n_cp)
    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        basisfunctions(i, u[j], degree, U, N[j, (i - degree):(degree + i)])

    b = np.copy(D[0:dim-1, 1:(n_cp-1)])

    N = np.delete(N, 0, 0)
    N = np.delete(N, -1, 0)
    N = np.delete(N, 0, 1)
    N = np.delete(N, -1, 1)

    dx = np.copy(D[0, :])
    dy = np.copy(D[1, :])

    b = N.T @ b.T
    N = N.T @ N
    p = la.lstsq(N, b)[0].T

    px = p[0, :]
    py = p[1, :]
    # insert bc
    px = np.insert(px, 0, bc0[0])
    px = np.append(px, bc1[0])
    py = np.insert(py, 0, bc0[1])
    py = np.append(py, bc1[1])

    # control points
    P[0, :] = px
    P[1, :] = py

    # points
    pts = 1000
    # curve points
    C = np.zeros((dim, pts))
    # basis functions
    N = np.zeros((pts, n_splines))

    u = np.linspace(min(U), max(U), pts)
    for j in range(pts):
        i = findspan(n_splines - 1, degree, u[j], U)
        nurbs(n_splines, degree, U, P, u[j], C[:, j], N[j, (i - degree):(degree + i)])

    cx = C[0, :]
    cy = C[1, :]

    dx = np.insert(dx, 4, 0.3)
    dy = np.insert(dy, 4, 1.0)
    dx = np.insert(dx, 5, 0.85)
    dy = np.insert(dy, 5, 0.95)

    plt.figure(1, figsize=(12, 10))
    plt.plot(dx, dy, 'o', ms=20)
    plt.plot(px, py, 'r--', linewidth=3)
    plt.plot(dx, dy, linewidth=3)
    plt.plot(cx, cy, linewidth=3)
    plt.ylim(min(dy) - 0.1, max(dy) + 0.1)
    plt.xlim(min(dx) - 0.2, max(dx) + 0.2)
    plt.axis('off')
    plt.tight_layout()

#    i = 0
#    for xy in zip(px, py):
#        plt.annotate('$\mathbf{P_{%d}}$'%i, xy=xy, textcoords='offset points', fontsize=52).draggable()
#        i +=1

    plt.savefig('controlpoints-interpolation.eps', transparent=True, dpi=300)
    plt.show()

def test_interpolation_curve():
    from curve import nurbs
    import numpy.linalg as la
    from findspan import findspan
    from basisfunctions import basisfunctions
    from knot_multiplicity import knot_multiplicity

    # weight
    w = 1.0
    # number of control pts
    n_cp = 14
    # degree
    degree = 2
    # number of pts
    pts = n_cp - degree + 1
    # knot vector
    mink = 0.0
    maxk = 1.0
    U = np.linspace(mink, maxk, pts)
    # knot multiplicity
    U = knot_multiplicity(degree, U)
    # dimensions
    dim = 3
    # n_splines is the number of B-spline curves
    n_splines = len(U) - degree - 1
    # points
    pts = n_splines
    # control points
    P = np.zeros((dim, n_cp))
    P[-1, :] = 1.0
    # interpolation points
    D = np.zeros((dim, n_cp))
    # curve points
    C = np.zeros((dim, pts))
    # basis functions
    N = np.zeros((n_splines, n_splines))

    assert(len(U) == (n_cp + degree + 1))

    # initialize P_i
    D[0:dim-1, 0] = np.array([-0.5, -0.5])
    D[0:dim-1, 1] = np.array([-2.0, 2.0])
    D[0:dim-1, 2] = np.array([-1.5, 3.0])
    D[0:dim-1, 3] = np.array([-0.5, 3.5])
    D[0:dim-1, 4] = np.array([0.0, 2.5])
    D[0:dim-1, 5] = np.array([-0.2, 1.0])
    D[0:dim-1, 6] = np.array([1.0, 0.95])
    D[0:dim-1, 7] = np.array([0.5, 3.0])
    D[0:dim-1, 8] = np.array([0.8, 5.0])
    D[0:dim-1, 9] = np.array([2.8, 5.2])
    D[0:dim-1, 10] = np.array([3.0, 1.2])
    D[0:dim-1, 11] = np.array([2.6, 0.0])
    D[0:dim-1, 12] = np.array([1.8, -0.8])
    D[0:dim-1, 13] = np.array([0.2, -0.2])

    u = np.linspace(mink, maxk, n_splines)
    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        basisfunctions(i, u[j], degree, U, N[j, (i - degree):(degree + i)])

    dx = np.copy(D[0, :])
    dy = np.copy(D[1, :])

    P[0, :] = la.solve(N, D[0, :])
    P[1, :] = la.solve(N, D[1, :])

    px = np.copy(P[0, :])
    py = np.copy(P[1, :])

    # points
    pts = 1000
    # curve points
    C = np.zeros((dim, pts))
    # basis functions
    N = np.zeros((pts, n_splines))

    u = np.linspace(min(U), max(U), pts)
    for j in range(pts):
        i = findspan(n_splines - 1, degree, u[j], U)
        nurbs(n_splines, degree, U, P, u[j], C[:, j], N[j, (i - degree):(degree + i)])

    cx = C[0, :]
    cy = C[1, :]

    plt.figure(1, figsize=(20, 14))
    plt.plot(u, N, linewidth=3)
    plt.plot(U, np.zeros_like(U), 'o', ms=25)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.tight_layout()

    plt.figure(2, figsize=(16, 14))
    plt.plot(dx, dy, 'o', ms=20)
    plt.plot(px, py, 'r--', linewidth=3)
    plt.plot(dx, dy, linewidth=3)
    plt.plot(cx, cy, linewidth=3)
    plt.xlim(1.1 * min(px), 1.1 * max(px))
    plt.ylim(1.1 * min(py), 1.1 * max(py))
    plt.axis('off')
    plt.tight_layout()

    i = 0
    for xy in zip(px, py):
        plt.annotate('$\mathbf{P_{%d}}$'%i, xy=xy, textcoords='offset points', fontsize=52).draggable()
        i +=1

    plt.savefig('controlpoints-interpolation.eps', transparent=True, dpi=300)
    plt.show()


# Michigan Tech Lecture Notes for CS 3621

def test_weighted_curve():
    from curve import nurbs
    from findspan import findspan
    from knot_multiplicity import knot_multiplicity

    # weight
    w = 2.0
    # number of control pts
    n_cp = 14
    # degree
    degree = 6
    # number of pts
    pts = n_cp - degree + 1
    # knot vector
    mink = 0.0
    maxk = 1.0
    U = np.linspace(mink, maxk, pts)
    # knot multiplicity
    U = knot_multiplicity(degree, U)
    # points
    pts = 1000
    # dimensions
    dim = 3
    # n_splines is the number of B-spline curves
    n_splines = len(U) - degree - 1
    # control points
    P = np.zeros((dim, n_cp))
    P[-1, :] = 1.0
    P[-1, 9] = w
    # curve points
    C = np.zeros((dim, pts))
    # array for basis functions
#    N = np.zeros(degree + 1)
    # basis functions
    N = np.zeros((pts, n_splines))

    assert(len(U) == (n_cp + degree + 1))

    # initialize P_i
    P[0:dim-1, 0] = np.array([-0.5, -0.5])
    P[0:dim-1, 1] = np.array([-2.0, 2.0])
    P[0:dim-1, 2] = np.array([-1.5, 3.0])
    P[0:dim-1, 3] = np.array([-0.5, 3.5])
    P[0:dim-1, 4] = np.array([0.0, 2.5])
    P[0:dim-1, 5] = np.array([-0.2, 1.0])
    P[0:dim-1, 6] = np.array([1.0, 0.95])
    P[0:dim-1, 7] = np.array([0.5, 3.0])
    P[0:dim-1, 8] = np.array([0.8, 5.0])
    P[0:dim-1, 9] = np.array([2.8, 5.2])
    P[0:dim-1, 10] = np.array([3.0, 1.2])
    P[0:dim-1, 11] = np.array([2.6, 0.0])
    P[0:dim-1, 12] = np.array([1.8, -0.8])
    P[0:dim-1, 13] = np.array([0.2, -0.2])

    px = np.copy(P[0, :])
    py = np.copy(P[1, :])

    P[0:dim-1, 9] = np.array([2.8, 5.2]) * w

    # random u value whom we are trying to find the ith interval it lies in
    u = np.linspace(min(U), max(U), pts)

    for j in range(pts):
        i = findspan(n_splines - 1, degree, u[j], U)
        nurbs(n_splines, degree, U, P, u[j], C[:, j], N[j, (i - degree):(degree + i)])

    cx = C[0, :]
    cy = C[1, :]

    plt.figure(1, figsize=(20, 14))
    plt.plot(u, N, linewidth=3)
    plt.plot(U, np.zeros_like(U), 'o', ms=25)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.tight_layout()
    if w < 1:
        plt.savefig('basisfunction-weight-0%d.eps'%(10*w), transparent=True, dpi=300)
    else:
        plt.savefig('basisfunction-weight-%d.eps'%w, transparent=True, dpi=300)

    plt.figure(2, figsize=(16, 14))
    plt.plot(px, py, 'o', ms=20)
    plt.plot(px[9], py[9], 'ro', ms=17)
    plt.plot(px, py, linewidth=3)
    plt.plot(cx, cy, linewidth=3)
    plt.xlim(1.1 * min(px), 1.1 * max(px))
    plt.axis('off')
    plt.title('Weight at $\mathbf{P_9} = %0.1f$'%w, fontsize=60)
    plt.tight_layout()

    i = 0
    for xy in zip(px, py):
        plt.annotate('$\mathbf{P_{%d}}$'%i, xy=xy, textcoords='offset points', fontsize=52)
        i +=1
    if w < 1:
        plt.savefig('controlpoints-weight-0%d.eps'%(10*w), transparent=True, dpi=300)
    else:
        plt.savefig('controlpoints-weight-%d.eps'%w, transparent=True, dpi=300)
#    plt.show()



# Textbook B-spline curve example
#   - quadratic

def test_convexity():
    from findspan import findspan
    from curve import nurbs

    # knot vector
    U = np.array([0., 0., 0., 1.0/5.0, 2.0/5.0, 3.0/5.0, 4.0/5.0, 1., 1., 1.])
    # degree
    degree = 2
    # points
    pts = 1000
    # dimensions
    dim = 3
    # n_splines is the number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of control points
    n_cp = 7
    # control points
    P = np.zeros((dim, n_cp))
    P[-1, :] = 1.0
    # curve points
    C = np.zeros((dim, pts))
    # array for basis functions
    N = np.zeros((pts, n_splines))

    assert(len(U) == (n_cp + degree + 1))

    # initialize P_i
    P[0:dim-1, 0] = np.array([0.0, 0.0])
    P[0:dim-1, 1] = np.array([-0.7, 2.0])
    P[0:dim-1, 2] = np.array([3.0, 2.5])
    P[0:dim-1, 3] = np.array([1.75, -0.25])
    P[0:dim-1, 4] = np.array([0.5, -3.0])
    P[0:dim-1, 5] = np.array([4.0, -2.5])
    P[0:dim-1, 6] = np.array([3.5, -0.5])

    u = np.linspace(min(U), max(U), pts)
    for j in range(pts):
        i = findspan(n_splines - 1, degree, u[j], U)
        nurbs(n_splines, degree, U, P, u[j], C[:, j], N[j, (i - degree):(degree + i)])

    px = P[0, :]
    py = P[1, :]

    cx = C[0, :]
    cy = C[1, :]

    i = 0
    for xy in zip(px, py):
        plt.annotate('$\mathbf{P_{%d}}$'%i, xy=xy, textcoords='offset points', fontsize=40)
        i +=1

    idx0 = np.where(u > 0.4)[0][0]
    idx1 = np.where(u > 0.6)[0][0]

    cxx = np.array([C[0, idx0], C[0, idx1]])
    cyy = np.array([C[1, idx0], C[1, idx1]])
    i = 0
    for xy in zip(cxx, cyy):
        if i == 0:
            plt.annotate('$\mathbf{C}(2/5)$', xy=xy, textcoords='offset points', fontsize=40)
        elif i == 1:
            plt.annotate('$\mathbf{C}(3/5)$', xy=xy, textcoords='offset points', fontsize=40)
        i += 1

    plt.plot(C[0, idx0], C[1, idx0], 'mo', ms=11)
    plt.plot(C[0, idx1], C[1, idx1], 'mo', ms=11)

    plt.plot(px, py, 'o', ms=9)
    plt.plot(px, py, linewidth=3)
    plt.plot(cx, cy, linewidth=3)
    plt.axis('off')
    plt.savefig('convexity-test.eps', transparent=True, dpi=300)
    plt.show()

# Textbook B-spline curve example
#   - cubic Bezier curve

def test_curve():
    from curve import curvepoint

    # knot vector
    U = np.array([0., 0., 0., 0., 1., 1., 1., 1.])
    # degree
    degree = 3
    # points
    pts = 100
    # dimensions
    dim = 2
    # n_splines is the number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of control points
    n_cp = 4
    # control points
    P = np.zeros((dim, n_cp))
    # curve points
    C = np.zeros((dim, pts))
    # array for basis functions
    N = np.zeros((1, degree + 1))

    assert(len(U) == (n_cp + degree + 1))
    # initialize P_i
    P[:, 0] = np.array([0.0, 0.0])
    P[:, 1] = np.array([1.0, 3.0])
    P[:, 2] = np.array([4.0, 4.0])
    P[:, 3] = np.array([5.0, 0.0])
    # random u value whom we are trying to find the ith interval it lies in
    u = np.linspace(min(U), max(U), pts)

    for j in range(pts):
        curvepoint(n_splines, degree, U, P, u[j], C[:, j], N)

    px = P[0, :]
    py = P[1, :]

    cx = C[0, :]
    cy = C[1, :]

    plt.plot(px, py, 'o', ms=9)
    plt.plot(px, py, linewidth=3)
    plt.plot(cx, cy, linewidth=3)
    plt.show()


def test_findspan():
    from findspan import findspan

    # knot vector
    U = np.array([0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.])
    # degree
    p = 2
    # n is the number of B-spline curves
    n = len(U) - p - 1
    # random u value whom we are trying to find the ith interval it lies in
    u = 2.5

    print('Number of B-spline curves: %d'%n)

    i = findspan(n - 1, p, u, U)
    print('ith interval: %d'%i)


def test_collocation_points():
    from findspan import findspan
    from derivatives import derivatives
    from collocation import collocation_points
    from knot_multiplicity import knot_multiplicity

    # degree
    degree = 2
    # knot vector
    mink = 0.0
    maxk = 10.0
    U = np.linspace(mink, maxk, 10)
    # knot multiplicity
    U = knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = degree
    # random u value whom we are trying to find the ith interval it lies in
    pts = 1000
    u = np.linspace(min(U), max(U), pts)
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((pts, n + 1, n_splines))
    # collocation points
    co = np.zeros((2, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    collocation_points(co, U, degree)

    plt.figure(1, figsize=(11, 8))
    plt.plot(u, ders[:, 0, :], linewidth=3)
    plt.plot(co[0, :], np.zeros_like(co[0, :]), 'o', ms=9)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()


# This routine checks if the discrete form of the product rule holds true.
# d(u(x) * v(x)) / dx ?= u(x) * dv(x) / dx + v(x) * d(u(x)) / dx

def test_product_rule():
    import numpy.linalg as la
    from findspan import findspan
    from derivatives import derivatives
    from collocation import collocation_points
    from knot_multiplicity import knot_multiplicity

    degree = 5
    pts = 50

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
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
    n = 1

    collocation_points(co, U, degree)

    u = co[0, :]
    ders = np.zeros((pts, n + 1, n_splines))
    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    B = ders[:, 0, :]
    D = ders[:, 1, :]

    x_velocity = np.sin(u)
    y_velocity = np.cos(u)
    xy_velocity = np.sin(u) * np.cos(u)
    sol = np.cos(2.0 * u)

    alpha_x = la.solve(B, x_velocity)
    alpha_y = la.solve(B, y_velocity)
    alpha_xy = la.solve(B, xy_velocity)

    expanded = x_velocity * (D @ alpha_y) + y_velocity * (D @ alpha_x)
    conservative = D @ alpha_xy

    Linfe = la.norm(sol - expanded, np.inf)
    Linfc = la.norm(sol - conservative, np.inf)

    plt.figure(1, figsize=(16, 12))
    plt.plot(u, sol, linewidth=3, label='solution')
    plt.plot(u, expanded, 'o', ms=9, label='expanded error = %e'%Linfe)
    plt.plot(u, conservative, 'x', ms=9, mew=3, label='conservative error = %e'%Linfc)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=36, loc='lower right')
    plt.show()


def test_basisfunctions(type_of_basis):
    from findspan import findspan
    from bspline import de_Boor
    from basisfunctions import basisfunctions
    from derivatives import derivatives

    if type_of_basis == 'non-periodic':
        # knot vector
        U = np.array([0., 0., 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1., 1.])
        # degree
        degree = 2
        # number of B-spline curves
        n_splines = len(U) - degree - 1
        # random u value whom we are trying to find the ith interval it lies in
        pts = 1000
        u = np.linspace(min(U), max(U), pts)
        # basis functions
        N = np.zeros((pts, n_splines))

        for j in range(len(u)):
            i = findspan(n_splines - 1, degree, u[j], U)
            basisfunctions(i, u[j], degree, U, N[j, (i - degree):(degree + i)])

        plt.rc('text', usetex=True)
        plt.figure(1, figsize=(20, 12))
        plt.plot(U, np.zeros_like(U), '|', ms=10)
        plt.plot(u, N, linewidth=5)
        plt.ylim(0.0, 1.05)
        plt.xlim(-0.02, 1.02)
        plt.xticks(fontsize=50)
        plt.yticks(fontsize=50)
        plt.savefig('collocation-pts.eps', transparent=True, dpi=300)
        plt.show()
    elif type_of_basis == 'periodic':
        # knot vector
        U = np.array([-3, -2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10, 11, 12., 13.])
        # degree
        degree = 3
        # number of B-spline curves
        n_splines = len(U) - degree - 1
        # random u value whom we are trying to find the ith interval it lies in
        pts = 1000
        u = np.linspace(min(U), max(U), pts)
        # basis functions
        N = np.zeros((pts, pts, degree+1))

        de_Boor(N, u, U, degree)

        plt.figure(1, figsize=(12, 8))
        plt.plot(U, np.zeros_like(U), '|', ms=10)
        plt.plot(u, N[:, :, degree], linewidth=3)
        plt.ylim(0.0, 1.2)
        plt.show()
    else:
        print('Erroneous input')
        quit()


def test_tanh_basisfunctions():
    from findspan import findspan
    from bspline import de_Boor
    from basisfunctions import basisfunctions
    from derivatives import derivatives

    # knot vector
    U = np.array([0., 0., 0., 0., 0., 0., 1., 2., 3., 4., 5.,\
            6., 7., 8., 9., 10, 10, 10., 10., 10., 10.])
    U = U/5.0
    # degree
    degree = 5
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # random u value whom we are trying to find the ith interval it lies in
    pts = 1000
    u = np.linspace(min(U), max(U), pts)
    # basis functions
    N = np.zeros((pts, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        basisfunctions(i, u[j], degree, U, N[j, (i - degree):(degree + i)])

    pts = 15
    x = np.linspace(0, pts, pts)
    y = np.linspace(0, pts, pts)
    co = np.tanh(2.0 * (2.0*(x+1-1)/(pts+1-1.0)-1.0))/np.tanh(2.0) + 1.0

    print(co)
    x= np.zeros((1, degree+1))
    a = np.array([0, 1, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 4, 5])

    for j in range(len(co)):
        i = findspan(pts - 1, degree, co[j], U)
        derivatives(i, co[j], degree, 0, U, x)
        y[j] = x[0, a[j]]
        print(x)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(20, 12))

    plt.plot(U, np.zeros_like(U), '|', ms=10)
    plt.plot(u, N, linewidth=5)
    plt.plot(co, y, 'ko', ms=13,linewidth=3)
    plt.ylim(0.0, 1.05)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    #plt.title('Degree 5 B-spline', fontsize=40)
    plt.savefig('tanh-collocation-pts.eps', transparent=True, dpi=300)
    plt.show()

def test_hyperbolic_tangent_mesh(pts):
    x = np.linspace(0, pts, pts)
    y = np.tanh(2.0 * (2.0*(x+1-1)/(pts+1-1.0)-1.0))/np.tanh(2.0) + 1.0

    plt.figure(1, figsize=(20, 12))
    plt.plot(x, y, linewidth=5)
    plt.ylim(-0.05, 2.05)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('Number of points', fontsize=50, labelpad=10)
    plt.ylabel('Coordinates', fontsize=50, labelpad=10)
    plt.title('Tanh collocation points', fontsize=50)
    plt.tight_layout()
    plt.savefig('tanh.eps', transparent=True, dpi=300)
    plt.show()


def test_derivatives():
    from findspan import findspan
    from derivatives import derivatives

    # knot vector
    U = np.array([0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.])
    U = np.array([0., 0., 0., 0.,  2., 4., 6, 8., 8., 8., 8.])
#    U = np.array([0., 0., 0., 0.,  2., 4., 4., 6., 6., 6., 8., 8., 8., 8.])
    # degree
    degree = 3
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = degree
    # random u value whom we are trying to find the ith interval it lies in
    pts = 1000
    u = np.linspace(min(U), max(U), pts)
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((pts, n + 1, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    plt.figure(1, figsize=(11, 8))
    U = np.array([-0.1, -0.05, 0.05, 0.1,  2., 4., 6, 7.9, 7.95, 8.05, 8.1])
    plt.plot(U, np.zeros_like(U), '|', ms=19, mew=2)
    plt.plot(u, ders[:, 1, :], linewidth=3)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.ylim(-0.2, 1.5)
    plt.xlim(-0.2, 8.2)
#    plt.savefig('original.svg', dpi=300)
    plt.show()


def test_recursion_derivatives():
    from findspan import findspan
    from recursion_derivatives import recursion_derivatives

    # knot vector
    U = np.array([0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.])
    U = np.array([0., 0., 0., 0.,  2., 4., 6, 8., 8., 8., 8.])
#    U = np.array([0., 0., 0., 0.,  2., 4., 4., 6., 6., 6., 8., 8., 8., 8.])
    # degree
    degree = 3
    # m is the highest index in U (m + 1) knots
    m = len(U) - 1
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = 0
    # random u value whom we are trying to find the ith interval it lies in
    pts = 1000
    u = np.linspace(min(U), max(U), pts)
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((pts, n + 1, n_splines))

    for j in range(len(u)):
        for i in range(n_splines):
            recursion_derivatives(degree, m , U, i, u[j], n, ders[j, :, i])

    print(ders[:, 0, :])
    plt.figure(1, figsize=(11, 8))
    U = np.array([-0.1, -0.05, 0.05, 0.1,  2., 4., 6, 7.9, 7.95, 8.05, 8.1])
    plt.plot(U, np.zeros_like(U), '|', ms=19, mew=2)
    plt.plot(u, ders[:, 0, :], linewidth=3)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.ylim(-0.2, 1.5)
    plt.xlim(-0.2, 8.2)
#    plt.savefig('original.svg', dpi=300)
    plt.show()


def test_compare_time(test):
    import csv

    if (test == 'degree'):
        p = np.arange(1, 21, 1)
        ptime_recursion = test_recursion_time(p, test)
        ptime = test_time(p, test)
        temp = zip(p, ptime, ptime_recursion)

    elif (test == 'knots'):
        knots = np.array([10, 100, 200, 500, 1000, 2000, 4000, 6000, 8000, 10000])
        ktime_recursion = test_recursion_time(knots, test)
        ktime = test_time(knots, test)
        temp = zip(knots, ktime, ktime_recursion)

    with open('output.csv', "w") as output:
        writer = csv.writer(output, lineterminator='\n', delimiter=',')
        for val in temp:
            writer.writerow(val)


def test_recursion_time(variant, test):
    import time
    from knot_multiplicity import knot_multiplicity
    from recursion_derivatives import recursion_derivatives

    # derivative
    n = 0

    if (test == 'degree'):
        p = variant
        t = np.zeros(len(p))

        # number of knots
        kpts = 1000
        for k in range(len(p)):
            t[k] = time.time()

            # knot vector
            U = np.linspace(0.0, 10.0, kpts)
            # multiplicity of knots according to degree
            U = knot_multiplicity(p[k], U)
            kpts = len(U)
            # m is the highest index in U (m + 1) knots
            m = len(U) - 1
            # number of B-spline curves
            n_splines = len(U) - p[k] - 1
            # number of points
            pts = kpts
            u = np.linspace(0.0, 10.0, pts)
            # basis functions
            N = np.zeros((pts, n + 1, n_splines))

            for j in range(len(u)):
                for i in range(n_splines):
                    recursion_derivatives(p[k], m , U, i, u[j], 0, N[j, :, i])

            t[k] = time.time() - t[k]

    elif (test == 'knots'):
        kpts = variant
        # degree
        p = 12

        t = np.zeros(len(kpts))
        for k in range(len(kpts)):
            t[k] = time.time()

            # knot vector
            U = np.linspace(0.0, 10.0, kpts[k])
            # multiplicity of knots according to degree
            U = knot_multiplicity(p, U)
            kpts[k] = len(U)
            # m is the highest index in U (m + 1) knots
            m = len(U) - 1
            # number of B-spline curves
            n_splines = len(U) - p - 1
            # number of points
            pts = kpts[k]
            u = np.linspace(0.0, 10.0, pts)
            # basis functions
            N = np.zeros((pts, n + 1, n_splines))

            for j in range(len(u)):
                for i in range(n_splines):
                    recursion_derivatives(p, m , U, i, u[j], 0, N[j, :, i])

            t[k] = time.time() - t[k]
            print('%d done'%kpts[k])

    return t


def test_time(variant, test):
    import time
    from findspan import findspan
    from basisfunctions import basisfunctions
    from knot_multiplicity import knot_multiplicity

    if (test == 'degree'):
        p = variant
        # number of knots
        kpts = 1000

        t = np.zeros(len(p))
        for k in range(len(p)):
            t[k] = time.time()

            # knot vector
            U = np.linspace(0.0, 10.0, kpts)
            # multiplicity of knots according to degree
            U = knot_multiplicity(p[k], U)
            kpts = len(U)
            # number of B-spline curves
            n_splines = len(U) - p[k] - 1
            # number of points
            pts = kpts
            u = np.linspace(0.0, 10.0, pts)
            # basis functions
            N = np.zeros((pts, n_splines))

            for j in range(len(u)):
                i = findspan(n_splines - 1, p[k], u[j], U)
                basisfunctions(i, u[j], p[k], U, N[j, (i - p[k]):(p[k] + i)])

            t[k] = time.time() - t[k]

    elif (test == 'knots'):
        kpts = variant
        # degree
        p = 12

        t = np.zeros(len(kpts))
        for k in range(len(kpts)):
            t[k] = time.time()

            # knot vector
            U = np.linspace(0.0, 10.0, kpts[k])
            # multiplicity of knots according to degree
            U = knot_multiplicity(p, U)
            kpts[k] = len(U)
            # number of B-spline curves
            n_splines = len(U) - p - 1
            # number of points
            pts = kpts[k]
            u = np.linspace(0.0, 10.0, pts)
            # basis functions
            N = np.zeros((pts, n_splines))

            for j in range(len(u)):
                i = findspan(n_splines - 1, p, u[j], U)
                basisfunctions(i, u[j], p, U, N[j, (i - p):(p + i)])

            t[k] = time.time() - t[k]
    return t


def test_eigenvalues():
    import numpy.linalg as la
    import scipy.linalg as sla
    from findspan import findspan
    from derivatives import derivatives
    from knot_multiplicity import knot_multiplicity

    # degree
    degree = 2
    # pts
    pts = 100
    # knot vector
    U = np.linspace(0.0, 1.0, pts)
    # multiplicity of knots according to degree
    U = knot_multiplicity(degree, U)
    # number of B-spline curves
    n_splines = len(U) - degree - 1
    # number of derivatives
    n = 2
    pts = n_splines
    u = np.linspace(min(U), max(U), pts)
    # der[i, j] is the kth derivative where
    # i is the number of distinct knots and j is the degree
    ders = np.zeros((pts, n + 1, n_splines))

    for j in range(len(u)):
        i = findspan(n_splines - 1, degree, u[j], U)
        derivatives(i, u[j], degree, n, U, ders[j, :, (i - degree):(degree + i)])

    eig, eigvec = sla.eig(ders[:, n, :])
    eig = np.sort(abs(eig))/max(np.abs(eig))
    x = np.arange(len(eig))/len(eig)


    plt.figure(1, figsize=(11, 8))
    plt.plot(x, eig, 'o')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()


def compact_loc(alpha, Erow, B8up, pts):
    import numpy.linalg as la

    a = 1/128.0 * (93.0 + 70 * alpha)
    b = 1.0/16.0 * (7.0 + 18.0 * alpha)
    c = 1.0/32.0 * (-7.0 + 14.0 *alpha)
    d = 1.0/16.0 * (1.0 - 2.0 * alpha)
    e = 1.0/128.0 * (-1.0 + 2.0 *alpha)

    Erow.fill(0.0); B8up.fill(0.0)

    Erow[-4] = e/2; Erow[-3] = d/2; Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a
    Erow[4] = e/2; Erow[3] = d/2; Erow[2] = c/2; Erow[1] = b/2

    a82 = np.array([1./2+alpha, (1./2+alpha)*2.0, 1./2+alpha, 0.0, 0.0, 0.0, 0.0])/2.
    a83 = np.array([-1/8.+alpha/4., 1./2+alpha, 2.0*(5./8+3*alpha/4.), 1./2+alpha, -1/8.+alpha/4., 0.0, 0.0])/2.
    a84 = np.array([1./32-alpha/16., -3./16+3*alpha/8., 15./32+17*alpha/16., 2.0*(11/16.+5.*alpha/8.), \
            15./32+17*alpha/16., -3./16+3*alpha/8., 1./32-alpha/16.])/2.

    for i in range(len(a84)):
        B8up[0, i] = a82[i]
        B8up[1, i] = a83[i]
        B8up[2, i] = a84[i]

    B8down = np.fliplr(np.flipud(B8up))
    explicit8loc = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit8 = np.zeros((pts, pts))

    for i in range(pts):
        explicit8loc[i, :] = np.roll(Erow, i)
        implicit8[i, :] = np.roll(Irow, i)

    explicit8loc[0, :]  = 0.0; explicit8loc[0, 0]   = 1.0
    explicit8loc[-1, :] = 0.0; explicit8loc[-1, -1] = 1.0
    implicit8[0, :]  = 0.0; implicit8[0, 0]   = 1.0
    implicit8[-1, :] = 0.0; implicit8[-1, -1] = 1.0
    explicit8loc[1:4, :] = B8up
    explicit8loc[-4:-1, :] = B8down

    A8loc = la.solve(implicit8, explicit8loc)
    eig = la.eigvals(A8loc)
    return eig


def test_eig_nonperiodic_compact(degree, pts):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 2

    kexactB = np.linspace(0.0, np.pi, pts-4)
    kexact = np.linspace(0.0, np.pi, pts)

    # B-spline
    fullu, ders, _ = nonperiodic_init(mink, maxk, pts, degree, deriv)
    B = ders[:, 0, :]
    I = np.eye(B.shape[0])
#    _ ,_ ,_ ,_ , P = \
#    smooth_spline_matrix_init(mink, maxk, degree, pts, 0.0, B, 2, 'nonperiodic', matrix_type='full')

    A = (B @ la.inv(B.T @ B + 1e-8*I) @ B.T)[1:-1, 1:-1]
    eig = la.eigvals(A)

    eig = np.sort(np.real(eig))[::-1]
    eig = (eig - np.min(eig))/np.ptp(eig)
    eig = eig[:-2]
    eig = (eig - np.min(eig))/np.ptp(eig)
    print('8th order B-spline, max=%.15f, min=%.15f'%(np.max(np.real(eig)), np.min(np.real(eig))))

    # 8th order compact, 8th boundary
    alpha=0.49
    a = 1/128.0 * (93.0 + 70 * alpha)
    b = 1.0/16.0 * (7.0 + 18.0 * alpha)
    c = 1.0/32.0 * (-7.0 + 14.0 *alpha)
    d = 1.0/16.0 * (1.0 - 2.0 * alpha)
    e = 1.0/128.0 * (-1.0 + 2.0 *alpha)

    Erow = np.zeros(pts)
    Erow[-4] = e/2; Erow[-3] = d/2; Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a
    Erow[4] = e/2; Erow[3] = d/2; Erow[2] = c/2; Erow[1] = b/2
    a82 = np.array([1./256 + 127*alpha/128., 31/32 + alpha/16., 7./64 + 25*alpha/32., -7/32. + 7*alpha/16.,\
            7*(5-10.*alpha)/128., -7/32. + 7.*alpha/16, 7/64.-7*alpha/32., -1/32.+alpha/16., 1/256.-alpha/128.])
    a83 = np.array([-1./256 + alpha/128., 1/32 + 15*alpha/16., 57/64. + 7.*alpha/32. ,7./32 + 9*alpha/16., \
            -7*(5-10.*alpha)/128., 7/32. - 7.*alpha/16, -7/64.+7*alpha/32, 1/32.-alpha/16., -1/256.+alpha/128.])
    a84 = np.array([1./256 - alpha/128., -1/32 + alpha/16., 7/64. + 25.*alpha/32. ,25./32 + 7*alpha/16., \
            35/128.+29.*alpha/64, -7/32. + 7.*alpha/16, 7/64.-7*alpha/32., -1/32.+alpha/16., 1/256.-alpha/128.])

    B8up = np.zeros((3, pts))
    for i in range(len(a82)):
        B8up[0, i] = a82[i]
        B8up[1, i] = a83[i]
        B8up[2, i] = a84[i]

    B8down = np.fliplr(np.flipud(B8up))
    explicit8 = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit8 = np.zeros((pts, pts))

    for i in range(pts):
        explicit8[i, :] = np.roll(Erow, i)
        implicit8[i, :] = np.roll(Irow, i)

    explicit8[0, :]  = 0.0; explicit8[0, 0]   = 1.0
    explicit8[-1, :] = 0.0; explicit8[-1, -1] = 1.0
    implicit8[0, :]  = 0.0; implicit8[0, 0]   = 1.0
    implicit8[-1, :] = 0.0; implicit8[-1, -1] = 1.0

    explicit8[:3, :] = B8up
    explicit8[-3:, :] = B8down

    A8 = la.solve(implicit8, explicit8)
    eig8 = la.eigvals(A8)
    imag = np.imag(eig8)[::-1]
    pos_imag_eig8 = imag[imag>=0]
    kexact_pos =kexact[np.where(imag>=0)]
    neg_imag_eig8 = imag[imag<=0]
    kexact_neg =kexact[np.where(imag<=0)]

    print('8th order compact, 8th boundary, maxR=%.15f, minR=%.15f, maxI=%.15f, minI=%.15f'\
    %(np.max(np.real(eig8)), np.min(np.real(eig8)), np.max(np.imag(eig8)), np.min(np.imag(eig8))))

    # 8th compact, 6th boundary
    a62 = np.array([1.0/64 + 31.*alpha/32., 29./32. + 3.*alpha/16., 15./64. + 17*alpha/32., \
        -5./16 + 5*alpha/32., 15./64 - 15.*alpha/32, -3./32 + 3*alpha/16., 1/64. - alpha/32.,0.0, 0.0])
    a63 = np.array([-1.0/64 + alpha/32., 3./32. + 13.*alpha/16., 49./64. + 15*alpha/32., \
        5./16 + 3*alpha/8., -15./64 + 15.*alpha/32, 3./32 - 3*alpha/16., -1/64. + alpha/32.,0.0, 0.0])
    a64 = a63.copy()

    B8up = np.zeros((3, pts))
    for i in range(len(a62)):
        B8up[0, i] = a62[i]
        B8up[1, i] = a63[i]
        B8up[2, i] = a64[i]

    B8down = np.fliplr(np.flipud(B8up))

    explicit8[0, :]  = 0.0; explicit8[0, 0]   = 1.0
    explicit8[-1, :] = 0.0; explicit8[-1, -1] = 1.0

    explicit8[:3, :] = B8up
    explicit8[-3:, :] = B8down

    A8 = la.solve(implicit8, explicit8)
    eig86 = la.eigvals(A8)

    print('8th order compact, 6th boundary, maxR=%.15f, minR=%.15f, maxI=%.15f, minI=%.15f'\
    %(np.max(np.real(eig86)), np.min(np.real(eig86)), np.max(np.imag(eig86)), np.min(np.imag(eig86))))

    # 8th compact, 4th boundary
    a42 = np.array([1.0/16.0 + 7*alpha/8., 3./4 + alpha/2., 3./8 + alpha/4., -1./4 + alpha/2., 1./16. - alpha/8., \
            0., 0., 0., 0.])
    a43 = a42.copy()
    a44 = a42.copy()

    B8up = np.zeros((3, pts))
    for i in range(len(a42)):
        B8up[0, i] = a42[i]
        B8up[1, i] = a43[i]
        B8up[2, i] = a44[i]

    B8down = np.fliplr(np.flipud(B8up))

    explicit8[0, :]  = 0.0; explicit8[0, 0]   = 1.0
    explicit8[-1, :] = 0.0; explicit8[-1, -1] = 1.0

    explicit8[:3, :] = B8up
    explicit8[-3:, :] = B8down

    A8 = la.solve(implicit8, explicit8)
    eig84 = la.eigvals(A8)

    print('8th order compact, 4th boundary, maxR=%.15f, minR=%.15f, maxI=%.15f, minI=%.15f'\
    %(np.max(np.real(eig84)), np.min(np.real(eig84)), np.max(np.imag(eig84)), np.min(np.imag(eig84))))

    # 8th order compact (low order centered)
    alpha=0.4
    eig8loc = compact_loc(alpha, Erow, B8up, pts)
    eig8loc = np.sort(np.real(eig8loc))[::-1]
    eig8loc = (eig8loc - np.min(eig8loc))/np.ptp(eig8loc)
    imagloc = np.imag(eig8loc)

    alpha=0.42
    eig84 = compact_loc(alpha, Erow, B8up, pts)
    eig84 = np.sort(np.real(eig84))[::-1]

    alpha=0.45
    eig86 = compact_loc(alpha, Erow, B8up, pts)
    eig86 = np.sort(np.real(eig86))[::-1]

    alpha=0.47
    eig8 = compact_loc(alpha, Erow, B8up, pts)
    eig8 = np.sort(np.real(eig8))[::-1]

    print('8th order compact loc, max=%.15f, min=%.15f'%(np.max(np.real(eig8loc)), np.min(np.real(eig8loc))))

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    color_list = plt.cm.Set3(np.linspace(0, 1, 5))

    # plot real
    plt.plot(kexactB, (np.real(eig)), 'b--',  linewidth=3,  label='8th-order B-spline ($B^{(6)}, \lambda^2 = 10^{-6}$)')
    plt.plot(kexact, np.sort((np.real(eig8)))[::-1], 'r--',  linewidth=3, label='8th-order boundary, $\\alpha=0.49$')
    plt.plot(kexact, np.sort((np.real(eig86)))[::-1], color=color_list[3], linestyle='--',  \
            linewidth=3, label='6th-order boundary, $\\alpha=0.49$')
    plt.plot(kexact, np.sort((np.real(eig84)))[::-1], color=color_list[1], linestyle='--',  \
            linewidth=3, label='4th-order boundary, $\\alpha=0.49$')
    plt.plot(kexact, (np.real(eig8loc)), 'g--',  linewidth=3, label='2nd-order centered, $\\alpha=0.49$')

    plt.plot(kexact, np.ones(len(kexact)), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=np.min(np.real(eig8))-0.1, top=1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('tf-real-compact-8-bspline.eps', transparent=True, dpi=300)
    plt.show()


def test_adv_diff_deriv_transfer_function(degree, pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from periodic_derivatives import periodic_init
    from wavenumber_periodic import periodic_fix
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 6

    # initialize B-spline operators
    u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)
    B = periodic_fix(ders[:, 0, :], degree)
    convection = periodic_fix(ders[:, 1, :], degree)
    diffusion = periodic_fix(ders[:, 2, :], degree)
    D8 = periodic_fix(ders[:, 6, :], degree)

    x = u[1:]
    kexact = x.copy()

#    # check B-spline operators
#    num = convection @ la.solve(B, np.sin(x))
#    print(la.norm(np.cos(x) - num, np.inf))
#    num = diffusion @ la.solve(B, np.sin(x))
#    print(la.norm(-np.sin(x) - num, np.inf))

    # 8th order compact
    pts = len(x)
    alpha = 0.49
    a = 1/128.0 * (93.0 + 70 * alpha)
    b = 1.0/16.0 * (7.0 + 18.0 * alpha)
    c = 1.0/32.0 * (-7.0 + 14.0 *alpha)
    d = 1.0/16.0 * (1.0 - 2.0 * alpha)
    e = 1.0/128.0 * (-1.0 + 2.0 *alpha)

    Erow = np.zeros(pts)
    Erow[-4] = e/2; Erow[-3] = d/2; Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a
    Erow[4] = e/2; Erow[3] = d/2; Erow[2] = c/2; Erow[1] = b/2
    explicit8 = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit8 = np.zeros((pts, pts))

    for i in range(pts):
        explicit8[i, :] = np.roll(Erow, i)
        implicit8[i, :] = np.roll(Irow, i)

#    # 8th derivative
#    Erow = np.zeros(pts)
#
#    Erow[-8] = 17311.0; Erow[-7] = -351616.0; Erow[-6] = 3434760.0; Erow[-5] = -21445760.0;
#    Erow[-4] = 95023460.0; Erow[-3] = -302537088.0; Erow[-2] = 689491768.0; Erow[-1] = -1126894720.0;
#    Erow[0] = 326523770.0
#    Erow[8] = 17311.0; Erow[7] = -351616.0; Erow[6] = 3434760.0; Erow[5] = -21445760.0;
#    Erow[4] = 95023460.0; Erow[3] = -302537088.0; Erow[2] = 689491768.0; Erow[1] = -1126894720.0;
#
#    #Erow[-4] = 1.0; Erow[-3] = -8.0; Erow[-2] = 28.0; Erow[-1] = -56.0; Erow[0] = 70.0
#    #Erow[4] = 1.0; Erow[3] = -8.0; Erow[2] = 28.0; Erow[1] = -56.0
#
#    #Erow[-1] = 1.0; Erow[0] = -2.0
#    #Erow[1] = 1.0
#    D8 = np.zeros((pts, pts))
#
#    for i in range(pts):
#        D8[i, :] = np.roll(Erow, i)
#
    dx = x[1] - x[0]
#    D8 = D8 / (3628800.0 * dx**8)
#
#    # check D8 operator
#    num = D8 @ np.sin(x)
#    plt.plot(x, num)
#    plt.show()
#    print(la.norm(-np.sin(x) - num, np.inf))
#    quit()

    transfer = np.zeros(pts)
    transfer2 = np.zeros(pts)
    transfer4 = np.zeros(pts)
    transfer6 = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    median=int(len(u)/2)
    lamlabel = ['10^{-6}']
    median = int(len(u)/2)
    for i in range(pts):
        original = np.exp(i* 1j* x)
        filtered = la.solve(implicit8, explicit8 @ original)

        input = -convection @ original + diffusion @ original
        output = -convection @ filtered + diffusion @ filtered
        output2 = -convection @ original + diffusion @ original + 0.1 * D8 @ original * dx**5
        output4 = -convection @ original + diffusion @ original + 0.1 * D8 @ original * dx**4.5
        output6 = -convection @ original + diffusion @ original + 0.1 * D8 @ original * dx**4.

        transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
        transfer2[i] = np.median(np.abs(input))/np.median(np.abs(output2))
        transfer4[i] = np.median(np.abs(input))/np.median(np.abs(output4))
        transfer6[i] = np.median(np.abs(input))/np.median(np.abs(output6))

    transfer[0] = 1.0; transfer2[0] = 1.0; transfer4[0] = 1.0; transfer6[0] = 1.0

    plt.plot(x[:median], transfer2[:median], linestyle='--', linewidth=3, label='Hyperviscosity 8th derivative, $C=0.1$')
    plt.plot(x[:median], transfer4[:median], linestyle='--', linewidth=3, label='Hyperviscosity 6th derivative, $C=0.1$')
    plt.plot(x[:median], transfer6[:median], linestyle='--', linewidth=3, label='Hyperviscosity 4th derivative, $C=0.1$')
    plt.plot(x[:median], transfer[:median], color='b', linewidth=3, label='8th order B-spline \
            ($B^{(6)}, \lambda^2 = %s)$'%lamlabel[0])

    plt.plot(kexact[:median], np.ones(median), 'k', linewidth=3, label='Reference')
    plt.ylim(bottom=-0.1, top=1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('Wavenumber', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.3)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    #plt.savefig('transfer-deriv-hyper-bspline.eps', transparent=True, dpi=300)
    plt.show()
    return


def test_adv_diff_C_transfer_function(degree, pts, periodic=True):
    import numpy.linalg as la
    import scipy.linalg as sla
    from smooth_spline import smooth_spline_matrix_init
    from periodic_derivatives import periodic_init
    from wavenumber_periodic import periodic_fix
    from collocation import nonperiodic_init

    # knot vector
    mink = 0.0
    maxk = 2.0 * np.pi
    # number of derivatives
    deriv = 6

    # initialize B-spline operators
    u, ders, U = periodic_init(mink, maxk, pts, degree, deriv)
    B = periodic_fix(ders[:, 0, :], degree)
    convection = periodic_fix(ders[:, 1, :], degree)
    diffusion = periodic_fix(ders[:, 2, :], degree)
    D8 = periodic_fix(ders[:, 6, :], degree)

    x = u[1:]
    kexact = x.copy()

#    # check B-spline operators
#    num = convection @ la.solve(B, np.sin(x))
#    print(la.norm(np.cos(x) - num, np.inf))
#    num = diffusion @ la.solve(B, np.sin(x))
#    print(la.norm(-np.sin(x) - num, np.inf))

    # 8th order compact
    pts = len(x)
    alpha = 0.49
    a = 1/128.0 * (93.0 + 70 * alpha)
    b = 1.0/16.0 * (7.0 + 18.0 * alpha)
    c = 1.0/32.0 * (-7.0 + 14.0 *alpha)
    d = 1.0/16.0 * (1.0 - 2.0 * alpha)
    e = 1.0/128.0 * (-1.0 + 2.0 *alpha)

    Erow = np.zeros(pts)
    Erow[-4] = e/2; Erow[-3] = d/2; Erow[-2] = c/2; Erow[-1] = b/2; Erow[0] = a
    Erow[4] = e/2; Erow[3] = d/2; Erow[2] = c/2; Erow[1] = b/2
    explicit8 = np.zeros((pts, pts))

    Irow = np.zeros(pts)
    Irow[-1] = alpha; Irow[0] = 1; Irow[1] = alpha
    implicit8 = np.zeros((pts, pts))

    for i in range(pts):
        explicit8[i, :] = np.roll(Erow, i)
        implicit8[i, :] = np.roll(Irow, i)

#    # 8th derivative
#    Erow = np.zeros(pts)
#
#    Erow[-8] = 17311.0; Erow[-7] = -351616.0; Erow[-6] = 3434760.0; Erow[-5] = -21445760.0;
#    Erow[-4] = 95023460.0; Erow[-3] = -302537088.0; Erow[-2] = 689491768.0; Erow[-1] = -1126894720.0;
#    Erow[0] = 326523770.0
#    Erow[8] = 17311.0; Erow[7] = -351616.0; Erow[6] = 3434760.0; Erow[5] = -21445760.0;
#    Erow[4] = 95023460.0; Erow[3] = -302537088.0; Erow[2] = 689491768.0; Erow[1] = -1126894720.0;
#
#    #Erow[-4] = 1.0; Erow[-3] = -8.0; Erow[-2] = 28.0; Erow[-1] = -56.0; Erow[0] = 70.0
#    #Erow[4] = 1.0; Erow[3] = -8.0; Erow[2] = 28.0; Erow[1] = -56.0
#
#    #Erow[-1] = 1.0; Erow[0] = -2.0
#    #Erow[1] = 1.0
#    D8 = np.zeros((pts, pts))
#
#    for i in range(pts):
#        D8[i, :] = np.roll(Erow, i)
#
    dx = x[1] - x[0]
#    D8 = D8 / (3628800.0 * dx**8)
#
#    # check D8 operator
#    num = D8 @ np.sin(x)
#    plt.plot(x, num)
#    plt.show()
#    print(la.norm(-np.sin(x) - num, np.inf))
#    quit()

    transfer = np.zeros(pts)
    transfer2 = np.zeros(pts)
    transfer4 = np.zeros(pts)
    transfer6 = np.zeros(pts)
    transfer8 = np.zeros(pts)
    transfer10 = np.zeros(pts)
    pade = np.zeros(pts)

    plt.rc('text', usetex=True)
    plt.figure(1, figsize=(16, 8))
    plt.rcParams["font.family"] = "Times New Roman"

    median=int(len(u)/2)
    lamlabel = ['10^{-6}']
    median = int(len(u)/2)
    for i in range(pts):
        original = np.exp(i* 1j* x)
        filtered = la.solve(implicit8, explicit8 @ original)

        input = -convection @ original + diffusion @ original
        output = -convection @ filtered + diffusion @ filtered
        output2 = -convection @ original + diffusion @ original + 1e-3 * D8 @ original * dx**4.5
        output4 = -convection @ original + diffusion @ original + 0.3 * D8 @ original * dx**4.5
        output6 = -convection @ original + diffusion @ original + 0.5 * D8 @ original * dx**4.5
        output8 = -convection @ original + diffusion @ original + 0.7 * D8 @ original * dx**4.5
        output10 = -convection @ original + diffusion @ original + 1.0 * D8 @ original * dx**4.5

        transfer[i] = np.median(np.abs(output))/np.median(np.abs(input))
        transfer2[i] = np.median(np.abs(input))/np.median(np.abs(output2))
        transfer4[i] = np.median(np.abs(input))/np.median(np.abs(output4))
        transfer6[i] = np.median(np.abs(input))/np.median(np.abs(output6))
        transfer8[i] = np.median(np.abs(input))/np.median(np.abs(output8))
        transfer10[i] = np.median(np.abs(input))/np.median(np.abs(output10))

    transfer[0] = 1.0; transfer2[0] = 1.0; transfer4[0] = 1.0; transfer6[0] = 1.0; transfer8[0] = 1.0; transfer10[0] = 1.0

    den = 3.0 + 2.0 * np.cos(kexact)
    for i in range(pts):
        num = 80.0  + 2.0 * (-60.0 * np.cos(kexact[i]) + 24. * np.cos(2 * kexact[i]) -4.0 * np.cos(3* kexact[i]))
        pade[i] = 1./256. * num / den[i]

    print(transfer2[median-1])
    transfer2[:median] = 1.0 - 1.0/(1.0-transfer2[median-1]) * (1-transfer2[:median])

    plt.plot(kexact[:median], 1-np.real(pade[:median]), 'r--', linewidth=3, label='Hyperviscosity compact 6th derivative')
    plt.plot(x[:median], transfer2[:median], 'b--', linewidth=3, label='Hyperviscosity B-spline 6th derivative')
    #plt.plot(x[:median], transfer4[:median], linestyle='--', linewidth=3, label='Hyperviscosity 6th derivative, $C=0.3$')
    #plt.plot(x[:median], transfer6[:median], linestyle='--', linewidth=3, label='Hyperviscosity 6th derivative, $C=0.5$')
    #plt.plot(x[:median], transfer8[:median], linestyle='--', linewidth=3, label='Hyperviscosity 6th derivative, $C=0.7$')
    #plt.plot(x[:median], transfer10[:median], linestyle='--', linewidth=3, label='Hyperviscosity 6th derivative, $C=1.0$')
    plt.plot(x[:median], transfer[:median], color='b', linewidth=3, label='B-spline filter with $B^{(6)}$ penalty matrix \
            and $\lambda^2 = %s$'%lamlabel[0])

    plt.plot(kexact[:median], np.ones(median), 'k', linewidth=3, label='Reference')
    #plt.ylim(bottom=-0.1, top=1.1)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('Transfer function (filtering)', labelpad=8, fontsize=30)
    plt.xlabel('$k\Delta x$', labelpad=8, fontsize=30)
    plt.gca().set_aspect(1.5)
    plt.legend(fontsize=28, framealpha=0.0)
    plt.tight_layout()
    plt.savefig('transfer-hyper-bspline.eps', transparent=True, dpi=300)
    plt.show()
    return


def main():
    #test_findspan()
    #test_basisfunctions('non-periodic')
    #test_hyperbolic_tangent_mesh(80)
    #test_tanh_basisfunctions()
    #test_derivatives()
    #test_time()
    #test_recursion_derivatives()
    #test_compare_time('knots')
    #test_eigenvalues()
    #test_collocation_points()
    #test_curve()
    #test_weighted_curve()
    #test_interpolation_curve()
    #test_approximation_curve()
    #test_convexity()
    #test_product_rule()
    #test_tikhonov_filter(5, 500)
    #test_smooth_spline_svd(5, 100)
    #test_smooth_spline_right_vector(5, 100)
    #test_smooth_spline_svd_coeff(5, 64)
    #test_smooth_spline_transfer_function_derivative(220)
    #test_smooth_spline_vs_compact_transfer(5, 100)
    #test_smooth_spline_transfer_function_lambda(220)
    #test_smooth_spline_vs_fd_transfer(5, 220)
    #test_smooth_spline_vs_fd_transfer_area(5, 220)
    #test_smooth_spline_vs_optimized_fd_transfer_area(5, 220)
    #test_boundary_fd_transfer(220)
    #test_eig_periodic_transfer(300)
    #test_eig_nonperiodic_transfer(9, 100)
    #test_derivative_order(7, 180)
    #output_derivative_order(7)
    #test_smooth_spline_vs_optimized_fd_transfer(5, 220)
    #test_smooth_spline_vs_drp7_optimized_transfer(5, 220)
    #test_ideal_transfer()
    #test_analytic_bspline_svd(5, 20)
    #plot_bspline_piecewise()
    #test_eig_nonperiodic_drp_9point_centered(5, 220)
    #test_eig_nonperiodic_drp_7point_antisym(5, 220)
    #test_eig_nonperiodic_drp(9, 220)
    #test_eig_nonperiodic_compact(9, 220)
    #test_channel_flow_ic()
    #test_adv_diff_deriv_transfer_function(7, 300)
    test_adv_diff_C_transfer_function(7, 300)


if __name__ == "__main__":
    main()
