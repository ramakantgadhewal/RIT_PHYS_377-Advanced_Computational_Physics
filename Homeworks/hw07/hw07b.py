"""
Author: Ramsey (Rayla) Phuc
"""
import matplotlib.pyplot as plt
import numpy as np


def pb01():
    """
    Part a: Replacing the second derivative in this equation with its finite-difference
    approximation (as given below), derive a relaxation-method equation for solving
    this problem on a time-like “grid” of points with separation h.
    """

    """
    x(t) = (1/2) * [x(t + h) + x(t - h) + g*h^2]
    """

    """
    Part b: Taking the boundary conditions to be that x = 0 at t = 0 and t = 10, 
    write a program to solve for the height of the ball as a function of time using 
    the relaxation method with 100 points and make a plot of the result from t = 0 to 
    t = 10. Run the relaxation method until the answers change by 10−6 or less at every 
    point on each step.
    """

    # Constants and Parameters
    g = 9.81

    # Initialize 2 arrays with size n
    n = 100
    x1, x2 = np.zeros(n, dtype=float), np.zeros(n, dtype=float)

    # Initializing the arrays
    x1[0] = x1[n - 1] = 0
    x1[1:n - 1] = 1
    x2[0] = x2[n - 1] = 0

    # Time interval
    t1, t2 = 0, 10
    ts = np.linspace(t1, t2, n)

    # Step size
    h = (t2 - t1) / n

    # Accuracy
    target = 10e-7

    while max(abs(x1 - x2)) > target:
        # Compute new position values
        for i in range(1, n - 1):
            x2[i] = (1 / 2) * (x1[i + 1] + x1[i - 1] + g * h ** 2)
        # Swap arrays
        x1, x2 = x2, x1
    plt.plot(ts, x1, label="x(t)")
    plt.legend()
    plt.savefig("hw07b.png")
    plt.show()
    return


def extraCredit():
    """

    """
    return


if __name__ == '__main__':
    pb01()
    pass
