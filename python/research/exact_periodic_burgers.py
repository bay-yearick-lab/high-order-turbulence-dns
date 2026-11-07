import numpy as np
import matplotlib.pyplot as plt


def coefficient_a(n, R):
    import scipy.special as ss
    return (-1)**n * ss.ive(n, 0.5 * R)


def exact_burgers(x, n, R, t):

    a = []
    for i in range(n):
        a.append(coefficient_a(i, R))

    num = 0.0
    for i in range(1, n):
        num += i * a[i] * np.exp(-(i**2 * t)) * np.sin(i * x)

    den = 0.0
    for i in range(1, n):
        den += a[i] * np.exp(-(i**2 * t)) * np.cos(i * x)


    return (4.0 * num) / (a[0] + 2.0 * den)


def main():

    pts = 32
    x = np.linspace(-np.pi, np.pi, pts)

    R = 20.0
    t = 0.075

    exact = exact_burgers(x, pts, R, t)

    plt.plot(x, exact)
    plt.show()


if __name__ == "__main__":
    main()
