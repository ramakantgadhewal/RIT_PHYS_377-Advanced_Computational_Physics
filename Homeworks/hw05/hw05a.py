"""
Author: Ramsey (Rayla) Phuc
File: hw0a.py
"""
import numpy as np


def bisection_steps(interval, p):
    """
    Determines the number of iterations of the Bisection Method required to
    get a root within p decimal places in an interval [a, b]
    """
    a, b = interval
    return np.ceil((np.log(2 * (b - a)) + p * np.log(10)) / np.log(2))


def binary(f, interval, n, target):
    a, b = interval
    if f(a) * f(b) > 0:
        raise Exception("Intermediate Value Theorem is not satisfied. "
                        "Bisection method will fail.")
    c = 0
    avals, bvals, xs = [a], [b], []
    while c <= n:
        x = (a + b) / 2
        xs.append(x)
        # print("x = " + str(x) + ", n = " + str(c))
        c += 1
        if abs(f(x)) < target:
            avals = np.asarray(avals, dtype=float)
            bvals = np.asarray(bvals, dtype=float)
            xs = np.asarray(xs, dtype=float)
            return avals, bvals, xs, c
        else:
            if f(a) * f(x) < 0:
                b = x
            else:
                a = x
            avals.append(a)
            bvals.append(b)
    raise Exception("Cannot find root using the Bisection Method")


def secant(f, interval, TOL):
    xs = interval
    x0, x1 = xs[0], xs[1]
    x = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0))
    c = 3
    while abs(f(x)) > TOL:
        x0 = x1
        x1 = x
        x = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0))
        c += 1
        xs.append(x)
    return np.asarray(xs, dtype=float), c


def pb01():
    def f(x):
        """
        Equation to find the root of
        """
        return 5 * np.exp(-x) + x - 5

    def I(wavelength):
        """
        Planck's radiation law
        """
        num = 2 * np.pi * h * (c ** 2) * (wavelength ** (-5))
        denom = np.exp(h * c / (wavelength * kb * T)) - 1

        return num / denom

    # Constants
    h = 6.62607004 * 10 ** (-34)  # Planck's Constant in m^2 kg / s
    c = 299792458  # Speed of Light in m/s
    kb = 1.38064852 * 10 ** (-23)  # Boltzmann constant
    wavelength = 502 * 10 ** (-9)  # nm -> m

    # Initial Interval
    x1, x2 = 1, 10
    interval = [x1, x2]

    # Target Accuracy
    p = 6
    target = 1 * 10 ** (-p)

    # Determine the number of steps needed
    n = bisection_steps(interval, p)

    # Find displacement constant
    x = binary(f, interval, n, target)[2][-1]
    print("Root of f(x) is: " + str(x))

    # Calculating Wien displacement constant
    b = h * c / (x * kb)
    print("Displacement Constant: " + str(b))

    # Calculating Surface temperature of sun
    # Wein's displacement law
    T = b / wavelength
    print("Surface temperature of the Sun is: " + str(T))
    return


def pb02():
    def f(r):
        """

        """
        return omega ** 2 * r - (G * M / (r ** 2)) + (G * m / ((R - r) ** 2))

    # Constants
    G = 6.674 * 10 ** (-11)
    M = 5.974 * 10 ** 24
    m = 7.348 * 10 ** 22
    R = 3.844 * 10 ** 8
    omega = 2.6622 * 10 ** (-6)

    # Initial Interval
    r1, r2 = 1, 10
    interval = [r1, r2]

    # Target Accuracy
    p = 4
    target = 1 * 10 ** (-p)

    # Find displacement constant
    x = secant(f, interval, target)[0][-1]
    print("Root of f(x) is: " + str(x))
    return


if __name__ == '__main__':
    print("================= Problem 01 =================")
    pb01()
    print("\n================= Problem 02 =================")
    pb02()
