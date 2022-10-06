"""
Author: Ramsey (Rayla) Phuc
"""

# Imported to use complex mathematical operations
import math

# Imported to create and display a plot
import matplotlib.pyplot as plt
# Imported to generate evenly spaced intervals
import numpy as np


def rk4_01(df, y, xs, RC):
    """
    RK4 from hw03a but with the following modifications:
        1. Adding an additional parameter RC
    """
    # Interval and number of steps
    a, b, n = xs[0], xs[-1], len(xs)
    # Step Size
    h = (b - a) / n
    # List of y values
    ys = [y]
    # For each x value
    for i in range(1, n):
        # Calculate k1 = h * df(y, x)
        k1 = h * df(ys[i - 1], xs[i], RC)
        # Calculate k2 = h * df(y + k1/2, x+ h/2)
        k2 = h * df(ys[i - 1] + k1 / 2, xs[i] + h / 2, RC)
        # Calculate k3 = h * df(y + k2/2, x + h/2)
        k3 = h * df(ys[i - 1] + k2 / 2, xs[i] + h / 2, RC)
        # Calculate k4 = h * df(y + k3, x + h)
        k4 = h * df(ys[i - 1] + k3, xs[i] + h, RC)
        # Calculate y(x+h) = y(x) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        y = ys[i - 1] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        # Add y(x + h) to ys
        ys.append(y)
    # Return ys as an array
    return np.asarray(ys)


def pb01():
    """
    Code for Problem 01
    """

    def V_in(t):
        if math.floor(2 * t) % 2 == 0:
            return 1
        else:
            return -1

    def df(V_out, t, RC):
        return (V_in(t) - V_out) / RC

    # Initial Conditions: f(t) = V_out
    V_out, t = 0, 0
    # time interval
    t1, t2 = 0, 10
    # Number of steps
    h = 0.001
    # Create time interval
    ts = np.linspace(t1, t2, int(1 / h) + 1)
    # Different values of RC
    RCs = [0.01, 0.1, 1]
    # Different colors
    colors = ['k', 'r', 'b']
    # For each RC value
    for i in range(0, len(RCs)):
        # Solve the differential equation
        sol = rk4_01(df, V_out, ts, RCs[i])
        # Plot the results
        plt.plot(ts, sol, color=colors[i], label="RC = " + str(RCs[i]))
    # Display the legend
    plt.legend()
    # Showing the data
    plt.show()


def runge_kutta_4(df, r, ts, params):
    """
    RK4 from hw03a but with the following modifications:
        1. Allowing the use of multiple parameters
        2. Being able to solve a system of differential equations instead of
           just one differential equation
    """
    # Interval and number of steps
    a, b, n = ts[0], ts[-1], len(ts)
    # Step Size
    h = (b - a) / n
    # List of y values
    rs = [r]
    # For each x value
    for i in range(1, n):
        # Calculate k1 = h * df(r, x)
        k1 = h * df(rs[i - 1], ts[i], params)
        # Calculate k2 = h * df(r + k1/2, x+ h/2)
        k2 = h * df(rs[i - 1] + k1 / 2, ts[i] + h / 2, params)
        # Calculate k3 = h * df(r + k2/2, x + h/2)
        k3 = h * df(rs[i - 1] + k2 / 2, ts[i] + h / 2, params)
        # Calculate k4 = h * df(r + k3, x + h)
        k4 = h * df(rs[i - 1] + k3, ts[i] + h, params)
        # Calculate y(r+h) = y(r) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        r = rs[i - 1] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        # Add y(x + h) to Fs
        rs.append(r)
    # Return rs as an array
    return np.asarray(rs)


def pb02():
    def df(F, t, params):
        """
        Creates a system of Differential Equations to solve numerically
        :param F: A tuple of Dependent Variables
        :param t: t value
        :param params: A tuple of parameters used in the differential equations
        :return: The system of Differential Equations
        """
        # Dependent Variables
        X, Y = F
        # Any defined parameters
        alpha, beta, gamma, delta = params
        # All the relevant differential equations
        dXdT = alpha * X - beta * X * Y
        dYdT = gamma * X * Y - delta * Y
        # Return the system of differential equations as an array
        return np.array([dXdT, dYdT])

    # Names for the corresponding independent variables
    vars = ['Prey', 'Predator']
    # Initial Conditions
    x0, y0 = 2, 2
    IC = np.array([x0, y0])
    # Parameters
    alpha, beta, gamma, delta = 1, 0.5, 0.5, 2
    params = alpha, beta, gamma, delta
    # Time Interval
    t1, t2 = 0, 30
    # Step Counter
    h = 0.01
    # Number of intervals
    n = int((t2 - t1) / h)
    # Create t values
    ts = np.linspace(t1, t2, n + 1)
    # Solve the differential equation
    sol = runge_kutta_4(df, IC, ts, params)

    # Different colors
    colors = ['k', 'r', 'b']

    for i in range(0, len(sol[0])):
        ys = sol[:, i]
        plt.plot(ts, ys, color=colors[i], label=vars[i])
    plt.legend()
    plt.show()
    """
    The graph represents the relationship between the Prey and Predator 
    populations. The oscillations shows how one population interacts with 
    the other. 
    """
    pass


if __name__ == '__main__':
    pb01()
    pb02()
    pass
