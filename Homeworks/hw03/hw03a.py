"""
Author: Ramsey (Rayla) Phuc
"""

# Imported to use complex mathematical operations
import math

# Imported to create and display a plot
import matplotlib.pyplot as plt
# Imported to generate evenly spaced intervals
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


def euler(f, x, t, ts):
    # Interval and number of steps
    a, b, n = ts[0], ts[-1], len(ts)

    # Step Size
    h = (b - a) / n

    # List of x values
    xs = [x]

    # For each t value
    for i in range(1, n):
        # Calculate x(t+h) + x(t) + h*f(x, t)
        x_n = xs[i - 1] + h * f(xs[i - 1], ts[i])

        # Add x(t+h) to xs
        xs.append(x_n)

    # Return xs as an array
    return np.asarray(xs)


def runge_kutta_2(f, x, t, ts):
    # Interval and number of steps
    a, b, n = ts[0], ts[-1], len(ts)

    # Step Size
    h = (b - a) / n

    # List of x values
    xs = [x]

    # For each t value
    for i in range(1, n):
        # Calculate k1 = h*f(x, t)
        k1 = h * f(xs[i - 1], ts[i])

        # Calculate k2 = h*f(x+k1/2, t+h/2)
        k2 = h * f(xs[i - 1] + k1 / 2, ts[i] + h / 2)

        # Calculate x(t+h) + x(t) + k2
        x_n = xs[i - 1] + k2

        # Add x(t+h) to xs
        xs.append(x_n)

    # Return xs as an array
    return np.asarray(xs)


def runge_kutta_4(f, x, t, ts):
    # Interval and number of steps
    a, b, n = ts[0], ts[-1], len(ts)

    # Step Size
    h = (b - a) / n

    # List of x values
    xs = [x]

    # For each t value
    for i in range(1, n):
        # Calculate k1 = h*f(x, t)
        k1 = h * f(xs[i - 1], ts[i])

        # Calculate k2 = h*f(x+k1/2, t+h/2)
        k2 = h * f(xs[i - 1] + k1 / 2, ts[i] + h / 2)

        # Calculate k3 = h*f(x+k2/2, t+h/2)
        k3 = h * f(xs[i - 1] + k2 / 2, ts[i] + h / 2)

        # Calculate k4 = h*f(x+k3, t+h)
        k4 = h * f(xs[i - 1] + k3, ts[i] + h)

        # Calculate x(t+h) = x(t) + (1/6)*(k1+2*k2+2*k3+k4)
        x_n = xs[i - 1] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

        # Add x(t+h) to xs
        xs.append(x_n)

    # Return xs as an array
    return np.asarray(xs)


def pb01():
    def f(x):
        return x ** 4 - 2 * x + 1

    # Interval and number of steps
    n1, n2 = 10, 20
    a, b = 0, 2

    # Integral results
    I1 = trapezoidal(f, a, b, n1)
    I2 = trapezoidal(f, a, b, n2)

    # Calculated Error via formula
    error = (1 / 3) * (I2 - I1)

    # Printing the results
    print("I = " + str(I2))
    print("Error via formula = " + str(error))
    print("Error via difference = " + str(4.4 - I2))

    """
    It is off by a factor c where in this case, c is approximately 2
    """

    pass


def pb02():
    def df(x, t):
        return -x ** 3 + math.sin(t)

    # Initial Conditions
    x, t = 0, 0

    # time interval
    t1, t2 = 0, 10

    # Number of steps
    n = 1000

    # Create time interval
    ts = np.linspace(t1, t2, n)

    # Solve the differential equation with the three methods
    solve1 = euler(df, x, t, ts)
    solve2 = runge_kutta_2(df, x, t, ts)
    solve3 = runge_kutta_4(df, x, t, ts)

    # Plot each method
    plt.plot(ts, solve1, 'k', label="Euler")
    plt.plot(ts, solve2, 'r', label="Runge-Kutta 2")
    plt.plot(ts, solve3, 'b', label="Runge-Kutta 4")

    # Add title to the plot
    plt.title("Computed using " + str(n) + " steps")

    # Display the legend
    plt.legend()

    # Showing the data
    print("Displaying plots")
    plt.show()

    """
    For smaller step sizes, we can see that there is a big difference between 
    the Runge-Kutta methods and Euler's method. However, there is a small 
    difference between the 2nd and 4th order Runge-Kutta methods.
    """
    pass


if __name__ == '__main__':
    print("-----Problem 1-----")
    # pb01()
    print("\n-----Problem 2-----")
    pb02()
    pass
