import numpy as np
import numpy.linalg as la
from smooth_spline import smooth_spline_matrix_init
from collocation import nonperiodic_init
import time

pts = 80
mink = 0.0
maxk = 2.0*np.pi
lamda =1.5e-6
degree = 5

u, ders, _ = nonperiodic_init(mink, maxk, pts, degree, 2)
B = ders[:, 0, :]
I = np.eye(B.shape[0])

_, _, _,_ , P = smooth_spline_matrix_init(mink, maxk, degree, pts, lamda, B, 2, 'nonperiodic')
lamda = np.random.rand(pts) * 1e-4
L = np.outer(lamda, lamda)
A = B.T@B + L * P

V1 = la.inv(B.T@B + 1e-8 * P)

# initialize operators
#B = np.random.rand(pts, pts)
#B = B.T @ B
#P = np.random.rand(pts, pts)
#P = P.T @ P
#lamda = np.random.rand(pts) * 1e-4
#L = np.outer(lamda, lamda)
#A = B + L * P
#I = np.eye(A.shape[0])

# find initial guess
#A1 = 0.0
#for i in range(pts):
#    A1 =  max(A1, np.sum(A[:, i]))
#Ainf = 0.0
#for i in range(pts):
#    Ainf =  max(Ainf, np.sum(A[i, :]))
#V1 = A.T / (A1 * Ainf)

tol = 1e-10
i = 0
total_time = 0.0
while la.norm(A@V1 - I, np.inf) > tol:
    start_time = time.time()
    AV1 = A @ V1
    V1 = V1 @ (2.0 * I - AV1)
    #V1 = 0.5 * V1 @ (9.0 * I - AV1 @ (16.0 * I -  AV1 @ (14.0 * I - AV1 @ (6*I - AV1))))
    total_time = total_time + time.time() - start_time
    i = i + 1
    print(i)

print(total_time, np.max(np.abs(la.inv(A) - V1)))
start_time = time.time()
Ainv = la.inv(B.T@B + np.diag(lamda)@ P @np.diag(lamda))
end_time = time.time() - start_time
print(end_time)
