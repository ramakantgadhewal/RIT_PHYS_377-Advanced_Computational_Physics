"""
Author: Ramsey (Rayla) Phuc
"""

import random

import numpy as np


def pb01():
    """
    Code for Problem 1
    """

    """
    Part a: Estimate the integral using the mean value method with 10000 points.
    """

    def f(x):
        """
        Function to integrate
        """
        return np.sin(1 / (x * (2 - x))) ** 2

    # Interval
    a, b = 0, 2

    # Number of points
    n = 10000

    # n evenly spaced points from a to b
    xs = np.linspace(a, b, n)

    # Remove the first and last elements
    xs = xs[1:-1]

    # Initialize the sums
    sum = 0

    # Perform Monte Carlo Integration (Mean Value Method)
    for i in range(0, n):
        sum += f(random.choice(xs))
    sum *= (b - a) / n

    # Print the result
    print("Approximated sum: " + str(sum))

    """
    Part b: Evaluate the error
    """
    # Initialize two sums
    sum2 = 0

    # Compute <f^2> and <f>^2
    for i in range(0, n):
        sum2 += (f(random.choice(xs))) ** 2
    sum = ((1 / n) * sum) ** 2
    sum2 = (1 / n) * sum2

    # Compute the variance
    variance = (b - a) * np.sqrt((sum2 - sum) / n)

    # Print the results
    print("Approximated Error: " + str(variance))
    return


def extraCredit():
    """
    Code for the extra credit problem
    """

    def f(r):
        """
        Returns 1 if the point is in the sphere. 0 otherwise
        """
        if np.sum(np.square(r)) <= 1:
            return 1
        else:
            return 0

    # v = np.array([1, 2, 3, 4, -5])
    # print(np.sum(np.square(v)))
    # print(np.dot(v, v))
    # interval of integration
    t1, t2 = -1, 1

    # Number of dimensions
    d = 10

    # Number of x points
    n = 1000000

    # List of n evenly spaced points from t1 to t2
    xs = np.linspace(t1, t2, n)

    # Initialize the sum
    sum = 0.0

    # Integrate over d-dimensions
    for i in range(0, n):
        r = np.zeros(d, dtype=float)
        for j in range(0, d):
            r[j] = random.choice(xs)
        sum += f(r)
    sum *= (2 ** d) / n

    # Print the result
    print(sum)


if __name__ == '__main__':
    pb01()
    extraCredit()
    pass
