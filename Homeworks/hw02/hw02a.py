"""
Author: Ramsey (Rayla) Phuc
File: hw02a.py
"""

# Imported to use complex mathematical operations
import math

# Imported to create and display a plot
import matplotlib.pyplot as plt
# Imported to read a .txt file
import numpy as np


def trapezoidal(f, a, b, n):
    """
    Performs the trapezoidal rule for approximating area
    :param f: f(x)
    :param a: lower bound
    :param b: upper bound
    :param n: Number of trapezoids
    :return: The area under f(x)
    """
    # Get the step interval
    h = (b - a) / n

    # Add half of the sum of f at the boundaries
    res = 0.5 * (f(a) + f(b))

    # Initialize i = 0
    i = 0

    # Add f(a + h * i) to the sum
    res += f(a + h * i)

    # Increment i by 1
    i += 1

    # While i < n
    while i < n:
        # Add f(a + h * i) to the sum
        res += f(a + h * i)
        # Increment i by 1
        i += 1
        pass

    # Multiply the sum by h
    res *= h

    # return the sum
    return res


def simpson(f, a, b, n):
    """
    Performs the simpson's rule for approximating area
    :param f: f(x)
    :param a: lower bound
    :param b: upper bound
    :param n: Number of trapezoids
    :return: The area under f(x)
    """
    # Get the interval step
    h = (b - a) / n

    # Add the value of f at the endpoints
    res = f(a) + f(b)

    # Doing summations for the odd case
    for i in range(1, n - 1, 2):
        res += 4 * f(a + i * h)
        pass

    # Doing summations for the even case
    for j in range(2, n - 2, 2):
        res += 2 * f(a + j * h)
        pass

    # multiply the current sum by h/3
    res *= h / 3

    # Return the result
    return res


def pb01():
    # Read the data
    data = np.loadtxt('velocities.txt')

    # Initialize t, x, and v
    t, x, v = [], [0], []

    # Extract the time and velocity values
    for e in data:
        t.append(e[0])
        v.append(e[1])
        pass

    # Set the sum = 0
    res = 0

    # Get the time interval
    a, b = int(t[0]), int(t[-1])

    # Compute each trapezoid and record the sum
    # Here, h = 1 so it is not needed when calculating
    for i in range(a, b):
        res += 0.5 * (v[i] + v[i + 1])
        x.append(res)
        pass

    # Print the sum
    print(res)

    # Plot both x and v vs t
    plt.plot(t, v, label="Velocity")
    plt.plot(t, x, label="Displacement")

    # Display the plot
    plt.legend()
    plt.show()
    pass


def pb02():
    def f(x):
        """
        Function to integrate
        """
        return x ** 4 - 2 * x + 1

    # Lower bound, Upper bound, and number of slices
    a, b, n = 0., 2., 10

    # Integrating via the Trapezoidal rule
    sum_t = trapezoidal(f, a, b, n)
    # Integrating via the Simpson's rule
    sum_s = simpson(f, a, b, n)

    # Print the results
    print("Area under f via Trapezoidal Rule: " + str(sum_t))
    print("Area under f via Simpson's Rule: " + str(sum_s))

    # Comparing the numerical results with the analytical result
    analytical_res = 4.4
    print("Fractional Error of Trapezoidal Rule: " + str(abs(analytical_res - sum_t)))
    print("Fractional Error of Simpson's Rule: " + str(abs(analytical_res - sum_s)))
    print("")

    # Integrating with different number of slices
    ns = [100, 1000]
    for n in ns:
        print("For n=" + str(n))
        # Integrating via the Trapezoidal rule
        sum_t = trapezoidal(f, a, b, n)

        # Print the results
        print("Area under f via Trapezoidal Rule: " + str(sum_t))
        print("")
        pass

    """
    The results gradually converges to the analytical answer with more steps.
    """

    pass


def pb03():
    def f(x):
        """
        Function to integrate
        """
        return math.exp(-x ** 2)

    # Initialize Lower bound, Upper bound, and step counter
    a, b, h = 0, 3, 0.1

    # Get the number of slices
    n = int((b - a) / h)

    # Initialize t and x
    t, x = [], [0]

    # Create values of t based on the initial values
    for i in range(0, n + 1):
        t.append(i * h)

    # Set the sum = 0
    res = 0

    # Compute each trapezoid and record the sum
    for i in range(0, n):
        res += (h / 2) * (f(a + (i - 1) * h) + f(a + i * h))
        x.append(res)
        pass

    # Print the sum
    print(res)

    # Plot x vs t
    plt.plot(t, x)
    plt.show()


if __name__ == '__main__':
    print("-----Problem 1-----")
    # pb01()
    print("\n-----Problem 2-----")
    pb02()
    print("\n-----Problem 3-----")
    # pb03()
    pass
