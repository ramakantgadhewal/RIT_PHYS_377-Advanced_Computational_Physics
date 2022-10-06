"""
Author: Ramsey (Rayla) Phuc
File: hw04b.py
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
    for i in range(0, len(ys)):
        plt.plot(xs, ys[i], color=COLORS[i], label=label[i])

    # Display the legend
    plt.legend()

    # Label the x and y axis
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    # Set the x and y axis boundaries. Uncomment to set boundaries.
    x_min, x_max = 0, 10
    y_min, y_max = -1, 1
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


def leapfrog(df, r, ts, h, params):
    # List of approximated solutions
    rs = [r]

    # List of midpoint values
    # Also calculate r(t + h/2) = r(t) + (h/2) * df(r, t) and append it to ms
    ms = [rs[0] + (h / 2) * df(rs[0], ts[0], params)]

    # For each t value
    for i in range(1, len(ts)):
        # Calculate r: r(t + h) = r(t) + h * df(r1, t + h/2)
        r = rs[i - 1] + h * df(ms[i - 1], ts[i - 1] + h / 2, params)

        # Calculate m: r(t + 3h/2) = r1 + h * df(r2, t + h)
        m = ms[i - 1] + h * df(r, ts[i - 1] + h, params)

        # Record the approximated r, m values
        rs.append(r)
        ms.append(m)
        pass

    return np.asarray(rs, dtype=float)


def verlet(df, F, ts, h, m, params):
    def gpe(params, m, r):
        """
        Gravitational Potential Energy
        """
        X, Y = r
        G, M = params
        return -G * M * m / np.sqrt(X ** 2 + Y ** 2)

    def ke(m, v):
        """
        Kinetic Energy
        """
        Vx, Vy = v
        return (1 / 2) * m * (Vx ** 2 + Vy ** 2)

    # List of approximated solutions
    rs, vs = [], []

    # Initialize the approximated solutions
    rs.append(F[0])
    vs.append(F[1])

    # List of midpoint values
    # Also calculate v(t + h/2) = v(t) + (h/2) * df(r, t) and append it to vs
    ms = vs[0] + (h / 2) * df([rs[0], vs[0]], ts[0], params)[-1]

    # List of Potential and Kinetic values
    PE, KE = [], []

    # Initialize the Potential and Kinetic Energies
    PE.append(gpe(params, m, rs[0]))
    KE.append(ke(m, vs[0]))

    # For each t value
    for i in range(1, len(ts)):
        # Calculate r: r(t + h) = r(t) + h * v(t + h/2)
        r = rs[i - 1] + h * ms

        # Calculate k = h * df(r(t + h), t + h)
        k = h * df([r, vs[i - 1]], ts[i - 1] + h, params)[-1]

        # Calculate v1: v(t + h) = v(t + h/2)  + (k/2)
        v = ms + (k / 2)

        # Calculate v2: v(t + 3h/2) = v(t + h/2)  + k
        ms += k

        # Calculate the Potential, Kinetic, and Total Energies
        E_potential = gpe(params, m, r)
        E_Kinetic = ke(m, v)

        # Record the approximated r, v values
        rs.append(r)
        vs.append(v)

        # Record the Potential, Kinetic, and Total Energies
        PE.append(E_potential)
        KE.append(E_Kinetic)

    # Convert each list to an array
    rs = np.asarray(rs, dtype=float)
    vs = np.asarray(vs, dtype=float)
    PE = np.asarray(PE, dtype=float)
    KE = np.asarray(KE, dtype=float)

    # Compute Total Energy
    TE = PE + KE

    return rs, vs, PE, KE, TE


def pb01():
    def df(r, t, params):
        # Dependent variables
        X, Y = r
        # Defined parameters
        g, l, C, omega = params
        # Differential equation(s) to solve
        dXdT = Y
        dYdT = Y ** 2 - X - 5
        return np.array([dXdT, dYdT], dtype=float)

    ############################################################################
    # Part a
    ############################################################################

    # Names for the corresponding independent variables
    variables = ['X']

    # Initial Conditions
    x0, y0 = 1, 0
    IC = np.array([x0, y0], dtype=float)

    # Parameters in meters, seconds
    g, l, C, omega = 9.81, 0.1, 2, 5
    params = g, l, C, omega

    # Time Interval
    t1, t2 = 0, 50

    # Step Size
    h = 0.001

    # Number of intervals
    n = int((t2 - t1) / h)

    # Create t values
    ts = np.linspace(t1, t2, n + 1)

    # Solve the differential equation
    sol = leapfrog(df, IC, ts, h, params)

    # Extract Data
    x_sol = sol[:, 0]
    xs = []
    for i in range(0, len(ts)):
        xs.append(x_sol[i])
    sols = [xs]
    # Plot data
    xaxis = "t"
    yaxis = "x(t)"
    title = ""
    plot_ode(ts, sols, variables, xaxis, yaxis, title)


def pb02():
    ############################################################################
    # Part (a)
    ############################################################################

    def df(F, t, params):
        # Dependent variables
        X, Y = F[0]
        Vx, Vy = F[1]
        # Defined parameters
        G, M = params
        r = np.sqrt(X ** 2 + Y ** 2)
        # Differential equation(s) to solve
        dXdT = Vx
        dVxdT = -G * M * X / (r ** 3)
        dYdT = Vy
        dVydT = -G * M * Y / (r ** 3)
        return np.array([[dXdT, dYdT], [dVxdT, dVydT]], dtype=float)

    ############################################################################
    # Part (b)
    ############################################################################

    # Initial Conditions
    distance = 1.471 * 10 ** 11  # m
    linear_velocity = 3.0287 * 10 ** 4  # m/s
    x0 = distance
    y0 = 0
    vx0 = 0
    vy0 = linear_velocity
    IC = np.array([
        [x0, y0],
        [vx0, vy0]
    ], dtype=float)
    # Parameters in meters, years
    G = 6.6738 * 10 ** (-11)  # m^3 kg^(-1) s^(-2)
    M = 1.9891 * 10 ** 30  # kg
    m = 5.9722 * 10 ** 24  # kg
    params = G, M

    # Time Interval
    t1, t2 = 0, 2 * (31557600)  # Years -> Seconds

    # Step Size
    h = 1 * (3600)  # Hours -> Seconds

    # Number of intervals
    n = int((t2 - t1) / h)

    # Create t values
    ts = np.linspace(t1, t2, n + 1)

    # Solve the differential equation
    sol = verlet(df, IC, ts, h, m, params)

    x_sol, y_sol = sol[0][:, 0], sol[0][:, 1]
    xs, ys = [], []
    for i in range(0, len(x_sol)):
        xs.append(x_sol[i])
        ys.append(y_sol[i])
    # Plot data
    xaxis = "time (s)"
    yaxis = "displacement (m)"
    title = ""
    variables = ['X', 'Y']
    sols = [ys]
    plot_ode(xs, sols, variables, xaxis, yaxis, title)

    ############################################################################
    # Part (c)
    ############################################################################

    PE, KE, TE = sol[2], sol[3], sol[4]
    yaxis = "Energy"
    variables = ["Potential", "Kinetic", "Total"]
    sols = [PE, KE, TE]
    plot_ode(ts, sols, variables, xaxis, yaxis, title)

    ############################################################################
    # Part (d)
    ############################################################################

    variables = ["Total"]
    sols = [TE]
    plot_ode(ts, sols, variables, xaxis, yaxis, title)
    pass


if __name__ == '__main__':
    # pb01()
    pb02()
