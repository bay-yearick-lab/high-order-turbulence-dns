import numpy as np

def knot_multiplicity(m, U):
    a = U[0]
    b = U[-1]

    # prepend/append
    for i in range(m):
       U = np.insert(U, 0, a)
       U = np.append(U, b)
    return U

def periodic_knot_multiplicity(m, U):
    dx = U[1] - U[0]

    # prepend/append
    for i in range(m):
       U = np.insert(U, 0, U[0] - dx)
       U = np.append(U, U[-1] + dx)
    return U
