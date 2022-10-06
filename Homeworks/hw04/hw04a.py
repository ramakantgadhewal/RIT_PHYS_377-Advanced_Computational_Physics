"""
Author: Ramsey (Rayla) Phuc
File: hw04a.py
"""
# Plotting
import matplotlib.pyplot as plt
# To create arrays and evenly spaced points
import numpy as np

# Names of certain colors and their HEX numbers
NAMES = ['Black', 'Red', 'Blue', 'Green', 'Cyan', 'Pink Lace',
         'Ice Blue', 'Light Hot Pink', 'Electric Blue', 'Deep Pink',
         'Fluorescent Blue', 'Hot Pink', 'Celeste', 'Rose Pink', 'Periwinkle',
         'Light Deep Pink', 'Ultra Pink']

COLORS = ['#000000', '#FF0000', '#0000FF', '#00FF00', '#00FFFF', '#FFDDF4',
          '#99FFFF', '#FFB3DE', '#7DF9FF', '#FF1493',
          '#15F4EE', '#FF69B4', '#B2FFFF', '#FF66CC', '#CCCCFF',
          '#FF5CCD', '#FF6FFF']


def plot_ode(xs, ys, label, xaxis, yaxis, title):
    """
    Generates a pretty graph based on the y values obtained from solving a
    system of ordinary differential equations
    :param xs: List of x values
    :param ys: A nested list of y values.
    :param label: Labels for each plot
    :param xaxis: Title for the x-axis
    :param yaxis: Title for the y-axis
    :param title: Title for the main graph
    """
    # Initialize a graph
    plt.figure()

    # Plot each data set
    if len(np.shape(ys)) > 1:
        for i in range(0, len(ys[0])):
            # Plot the data with a defined color and give the data a label
            plt.plot(xs, ys[:, i], color=COLORS[i], label=label[i])
    else:
        plt.plot(xs, ys, color=COLORS[0], label=label[0])

    # Display the legend
    plt.legend()

    # Label the x and y axis
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    # Set the x and y axis boundaries. Uncomment to set boundaries.
    x_min, x_max = 0, 10
    y_min, y_max = -2, 2
    # plt.xlim(x_min, x_max)
    # plt.ylim(y_min, y_max)

    # Give the graph a name
    plt.title(title)

    # Add grid lines to the graph
    plt.minorticks_on()
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')

    # Display the graph
    plt.show()
    pass


def rk4(df, r, ts, params):
    """
    Runge-Kutta 4 Method
    :param df: A system of differential equations to solve numerically
    :param r: A list of initial conditions. r = [x0, y0, z0, ...]
    :param ts: t values
    :param params: A tuple of parameters for the system of differential
                   equations
    :return:
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
        k1 = h * df(rs[i - 1], ts[i - 1], params)
        # Calculate k2 = h * df(r + k1/2, x+ h/2)
        k2 = h * df(rs[i - 1] + k1 / 2, ts[i - 1] + h / 2, params)
        # Calculate k3 = h * df(r + k2/2, x + h/2)
        k3 = h * df(rs[i - 1] + k2 / 2, ts[i - 1] + h / 2, params)
        # Calculate k4 = h * df(r + k3, x + h)
        k4 = h * df(rs[i - 1] + k3, ts[i - 1] + h, params)
        # Calculate y(r+h) = y(r) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        r = rs[i - 1] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        # Add y(r + h) to rs
        rs.append(r)
    # Return rs as an array
    return np.asarray(rs, dtype=object)


def pb01():
    def df(F, t, params):
        # Dependent variables
        theta, Y = F
        # Defined parameters
        g, l, C, omega = params
        # Differential equation(s) to solve
        dTHETAdT = Y
        dYdT = -(g / l) * np.sin(theta) + C * np.cos(theta) * np.sin(omega * t)
        return np.array([dTHETAdT, dYdT])

    ############################################################################
    # Part a
    ############################################################################

    # Names for the corresponding independent variables
    vars = ['theta', 'dTheta']
    # Initial Conditions
    theta, dTheta = 0, 0
    IC = np.array([theta, dTheta])
    # Parameters in meters, seconds
    g, l, C, omega = 9.81, 0.1, 2, 5
    params = g, l, C, omega
    # Time Interval
    t1, t2 = 0, 100
    # Number of intervals
    n = 10000
    # Create t values
    ts = np.linspace(t1, t2, n + 1)
    # Solve the differential equation
    sol = rk4(df, IC, ts, params)
    # Plot data
    xaxis = "x-axis"
    yaxis = "y-axis"
    title = ""
    plot_ode(ts, sol, vars, xaxis, yaxis, title)

    ############################################################################
    # Part B
    ############################################################################

    # Parameters
    omegas = np.linspace(omega * 0.9, omega * 1.1, 11)
    for omega in omegas:
        params = g, l, C, omega
        # Solve the differential equation
        sol = rk4(df, IC, ts, params)
        # Plot data
        xaxis = "x-axis"
        yaxis = "y-axis"
        title = "omega = " + str(omega)
        plot_ode(ts, sol, vars, xaxis, yaxis, title)
    pass


def pb02():
    def df1(F, t, params):
        X, Y = F
        omega = params
        dXdT = Y
        dYdT = -omega ** 2 * X
        return np.array([dXdT, dYdT])

    def df2(F, t, params):
        X, Y = F
        omega = params
        dXdT = Y
        dYdT = -omega ** 2 * X ** 3
        return np.array([dXdT, dYdT])

    ############################################################################
    # Part a
    ############################################################################

    # Names for the corresponding independent variables
    vars = ['x', 'dx']
    # Initial Conditions
    x, dx = 1, 0
    IC = np.array([x, dx])
    # Parameters in meters, seconds
    omega = 1
    params = omega
    # Time Interval
    t1, t2 = 0, 50
    # Number of intervals
    n = 1000
    # Create t values
    ts = np.linspace(t1, t2, n + 1)
    # Solve the differential equation
    sol = rk4(df1, IC, ts, params)
    # Plot data
    xaxis = "time (s)"
    yaxis = "Displacement (m)"
    title = "Problem 2, Part a"
    plot_ode(ts, sol, vars, xaxis, yaxis, title)

    ############################################################################
    # Part B: Increasing Amplitude
    ############################################################################

    # Initial Conditions
    x = 2
    IC = np.array([x, dx])
    # Solve the differential equation
    sol = rk4(df1, IC, ts, params)
    # Plot data
    xaxis = "time (s)"
    yaxis = "Displacement (m)"
    title = "Problem 2, Part b"
    plot_ode(ts, sol, vars, xaxis, yaxis, title)

    ############################################################################
    # Part C: Changing Differential Equation
    ############################################################################

    # Names for the corresponding independent variables
    vars = ['x', 'dx']
    # Initial Conditions
    x, dx = 1, 0
    IC = np.array([x, dx])
    # Parameters in meters, seconds
    omega = 1
    params = omega
    # Time Interval
    t1, t2 = 0, 50
    # Number of intervals
    n = 1000
    # Create t values
    ts = np.linspace(t1, t2, n + 1)
    # Solve the differential equation
    sol = rk4(df2, IC, ts, params)
    # Plot data
    xaxis = "time (s)"
    yaxis = "Displacement (m)"
    title = "Problem 2, Part c"
    plot_ode(ts, sol, vars, xaxis, yaxis, title)

    # Changing Amplitude
    # Initial Conditions
    x = 2
    IC = np.array([x, dx])
    # Solve the differential equation
    sol = rk4(df2, IC, ts, params)
    # Plot data
    xaxis = "time (s)"
    yaxis = "Displacement (m)"
    title = "Problem 2, Part c (changing amplitude)"
    plot_ode(ts, sol, vars, xaxis, yaxis, title)

    ############################################################################
    # Part D: Phase Space plot
    ############################################################################

    # Names for the corresponding independent variables
    vars = ['x', 'dx']
    # Initial Conditions
    x, dx = 1, 0
    IC = np.array([x, dx])
    # Parameters in meters, seconds
    omega = 1
    params = omega
    # Time Interval
    t1, t2 = 0, 50
    # Number of intervals
    n = 1000
    # Create t values
    ts = np.linspace(t1, t2, n + 1)
    # Solve the differential equation
    sol = rk4(df1, IC, ts, params)
    # Plot data
    xaxis = ""
    yaxis = ""
    title = "Problem 2, Part d"
    plot_ode(sol[:, 0], sol[:, 1], vars, xaxis, yaxis, "x vs dx")


def extra_credit():
    def df(r, t, params):
        X, Y = r
        mu, omega = params
        dXdT = Y
        dYdT = -omega ** 2 * X + mu * (1 - X ** 2) * Y
        return np.array([dXdT, dYdT])

    # Names for the corresponding independent variables
    vars = ['x', 'dx']
    # Initial Conditions
    x, dx = 1, 0
    IC = np.array([x, dx])
    # Parameters in meters, seconds
    mu, omega = 1, 1
    params = mu, omega
    # Time Interval
    t1, t2 = 0, 20
    # Number of intervals
    n = 10000
    # Create t values
    ts = np.linspace(t1, t2, n + 1)
    # Solve the differential equation
    sol = rk4(df, IC, ts, params)
    # Plot data
    xaxis = ""
    yaxis = ""
    title = "x vs dx"
    plot_ode(sol[:, 0], sol[:, 1], vars, xaxis, yaxis, title)

    # Changing mu value
    mus = [2, 4]
    for mu in mus:
        params = mu, omega
        # Solve the differential equation
        sol = rk4(df, IC, ts, params)
        title = "mu = " + str(mu)
        plot_ode(sol[:, 0], sol[:, 1], vars, xaxis, yaxis, title)
    pass


if __name__ == '__main__':
    pb01()
    pb02()
    extra_credit()
    pass
