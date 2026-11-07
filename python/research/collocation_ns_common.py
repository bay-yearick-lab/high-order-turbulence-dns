import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import matplotlib.pyplot as plt
from collocation import nonperiodic_init
from ns_mms import *


def low_storage_rk(stage):
    a = np.zeros((stage+1, stage+1))
    b = np.zeros(stage+1)
    c = np.zeros(stage+1)

    if (stage == 3):
        a[2, 1] = 8.0/15.0
        a[3, 1] = 1.0/4.0
        a[3, 2] = 5.0/12.0

        b[1] = 1.0/4.0
        b[3] = 3.0/4.0

        c[2] = 8.0/15.0
        c[3] = 2.0/3.0

    return a, b, c


def runge_kutta(stage):
    if (stage == 3):
        return np.array([1.0/3.0, 0.0, 2.0/3.0])
    elif (stage == 4):
        return np.array([0.0, 2.0/3.0, 0.0, 1.0/3.0])
    elif (stage == 6):
        return np.array([1.0/9.0, 2.0/5.0, 0.0, 4.0/9.0, 0.0, 2.0/45.0])
    elif (stage == 8):
        return np.array([2.0/15.0, 2.0/7.0, 2.0/9.0, 4.0/15.0, 0.0,
                        4.0/45.0, 0.0, 1.0/315.0])
    elif (stage == 10):
        return np.array([71./525, 22./81, 4./15, 4./21, 2./27, 4./75,
                        0.0, 8./945, 0.0, 2./14175])


def slice_matrix(matrix, row, column):
    index = np.ix_(row, column)
    return matrix[index]


def slice_matrix_1d(matrix):
    return matrix[1:-1, 1:-1]


def xboundary(f, t, fullu, xbc_index):

    nx = len(fullu[0])
    ny = len(fullu[1])

    x, y = np.unravel_index(xbc_index, (nx, ny))
    return f(fullu[0][x], fullu[1][y], t).flatten('F')


def yboundary(f, t, fullu, ybc_index):

    nx = len(fullu[0])
    ny = len(fullu[1])

    x, y = np.unravel_index(ybc_index, (nx, ny))
    return f(fullu[0][x], fullu[1][y], t).flatten('F')


def construct_full_solution(alpha, xbc, ybc, fullu, xbc_type):

    if xbc_type == 'dirichlet':
        nx = len(fullu[0])
        ny = len(fullu[1])

        b = np.append(xbc[0:ny], alpha)
        b = np.append(b, xbc[ny:])

        k = 0
        for i in range(nx-2):
            b = np.insert(b, (i+1)*ny    , ybc[k])
            b = np.insert(b, (i+2)*ny-1, ybc[k+1])
            k += 2

    elif xbc_type == 'periodic':
        nx = len(fullu[0])
        ny = len(fullu[1])

        b = alpha.copy()

        k = 0
        for i in range(nx):
            b = np.insert(b, (i)*ny    , ybc[k])
            b = np.insert(b, (i+1)*ny-1, ybc[k+1])
            k += 2

        b = np.append(b, b[0:ny])

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


def indexing(fullu, string, xbc_type):

    nx = len(fullu[0])
    ny = len(fullu[1])
    all_index = np.arange(nx*ny, dtype=int)

    if xbc_type == 'dirichlet':

        xbc_index = np.zeros(2*ny, dtype=int)
        ybc_index = np.zeros(2*(nx-2), dtype=int)

        # xbc row/column index
        xbc_index[:ny] = all_index[:ny]
        xbc_index[-ny:] = all_index[-ny:]

        # ybc row/column index
        k = 0
        for i in range(1, nx-1):
            ybc_index[k] = i * ny
            ybc_index[k+1] = (i + 1) * ny - 1
            k += 2

    elif xbc_type == 'periodic':

        xbc_index = None
        ybc_index = np.zeros(2*nx, dtype=int)

        # ybc row/column index
        k = 0
        for i in range(nx):
            ybc_index[k] = i * ny
            ybc_index[k+1] = (i + 1) * ny - 1
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
