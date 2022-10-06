"""
Author: Ramsey (Rayla) Phuc

Final Project Topic: Chaotic Systems - Double Pendulum
"""

"""
Importing custom modules
"""

"""
Importing built-in modules
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

"""
Global Variables
"""
NAMES = ["Black", "Red", "Blue", "Green", "Cyan", "Magenta",
         "Yellow", "Violet"]

COLORS = ["#000000", "#FF0000", "#0000FF", "#00FF00", "#00FFFF", "#FF00FF",
          "#FFFF00", "#7F00FF"]


################################################################################
# Helper functions
################################################################################
def plot(xs, ys, label, xaxis, yaxis, title, filename):
    """
    Generates a pretty graph based on the y values obtained

    :param xs: List of x values
    :param ys: A nested list of y values.
    :param label: List of labels
    :param xaxis: Title for the x-axis
    :param yaxis: Title for the y-axis
    :param title: Title for the main graph
    :param filename: save plot as <filename>.png
    """
    fig = plt.figure()
    # Plot each data set
    for i in range(0, len(ys)):
        if ys[i] is not None:
            if label[i] == "f":
                plt.plot(xs, ys[i], color=COLORS[i], label=label[i],
                         linewidth=1, markersize=1, marker="o")
            else:
                plt.plot(xs, ys[i], color=COLORS[i], label=label[i],
                         linewidth=1)
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
    plt.savefig(filename + '.png')
    plt.close()
    # plt.plot()
    pass


def get_data(data):
    """
    Transposes the array
    """
    numRows, numCols = len(data), len(data[0])
    temp = np.empty((numCols, numRows), dtype=float)
    for r in range(0, numRows):
        for c in range(0, numCols):
            temp[c][r] = data[r][c]
    if len(temp) == 1:
        return temp[0]
    else:
        return temp


def rk4(df, r, ts, h, params):
    """
    Runge-Kutta 4 Method

    :param df: System of ODE(s) to solve
    :param r: Initial Condition vector
    :param ts: time
    :param h: Step Size
    :param params: Parameters/Constants in the differential equation
    :return: Approximated solution
    """
    # Initialize the approximated solution
    rs = [r]

    # For each x value
    for i in range(1, len(ts)):
        # Compute f1 = f(r(t), t)
        f1 = df(rs[i - 1], ts[i - 1], params)

        # Compute f2 = f(r(t) + (h/2)*f1, t + (h/2))
        f2 = df(rs[i - 1] + (h / 2) * f1, ts[i - 1] + (h / 2), params)

        # Compute f3 = f(r(t) + (h/2)*f2, t + (h/2))
        f3 = df(rs[i - 1] + (h / 2) * f2, ts[i - 1] + (h / 2), params)

        # Compute f4 = f(r(t) + h*f3, t + h)
        f4 = df(rs[i - 1] + h * f3, ts[i - 1] + h, params)

        # Compute r(t + h) = r(t) + (h/6)*(f1 + 2*f2 + 2*f3 + f4)
        r = rs[i - 1] + (h / 6) * (f1 + 2 * f2 + 2 * f3 + f4)

        # Record the approximated value
        rs.append(r)

    # Return rs as an array
    return get_data(rs)


################################################################################
# Main function
################################################################################
def solve(m1, m2, l1, l2, g, theta1_0, theta2_0, omega1_0, omega2_0, t1, t2):
    """
    Solve the double pendulum system
    """

    def df(F, t, params):
        # Unpack variables
        theta1, theta2, omega1, omega2 = F

        # Unpack parameters
        m1, m2, l1, l2, g = params

        # System of differential equations
        theta1_t = omega1
        theta2_t = omega2

        # Breaking down each component since the expression is too big
        A1 = -g * (2 * m1 + m2) * np.sin(theta1)
        A2 = - m2 * g * np.sin(theta1 - 2 * theta2)
        A3 = - 2 * m2 * np.sin(theta1 - theta2)
        A4 = (omega2 ** 2) * l2 + (omega1 ** 2) * l1 * np.cos(theta1 - theta2)
        B = l1 * (2 * m1 + m2 - m2 * np.cos(2 * (theta1 - theta2)))
        omega1_t = (A1 + A2 + A3 * A4) / B

        # Breaking down each component since the expression is too big
        C1 = 2 * np.sin(theta1 - theta2)
        C2 = (omega1 ** 2) * l1 * (m1 + m2)
        C3 = g * (m1 + m2) * np.cos(theta1)
        C4 = (omega2 ** 2) * l2 * m2 * np.cos(theta1 - theta2)
        D = l2 * (2 * m1 + m2 - m2 * np.cos(2 * (theta1 - theta2)))
        omega2_t = C1 * (C2 + C3 + C4) / D
        return np.array([theta1_t, theta2_t, omega1_t, omega2_t], dtype=float)

    # Step size
    h = 0.0005

    # Initial Condition
    IC = np.array([np.radians(theta1_0), np.radians(theta2_0),
                   np.radians(omega1_0), np.radians(omega2_0)],
                  dtype=float)

    # Time span in seconds
    ts = np.arange(t1, t2, h)

    # Solve the model
    params = m1, m2, l1, l2, g
    theta1, theta2, omega1, omega2 = rk4(df, IC, ts, h, params)

    # Compute bob1's coordinates and velocity in Cartesian
    x1 = l1 * np.sin(theta1)
    y1 = - l1 * np.cos(theta1)
    vx1 = l1 * omega1 * np.cos(theta1)
    vy1 = l1 * omega1 * np.sin(theta1)
    v1 = vx1 ** 2 + vy1 ** 2

    # Compute bob2's coordinates and velocity in Cartesian
    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)
    vx2 = vx1 + l2 * omega2 * np.cos(theta2)
    vy2 = vy1 + l2 * omega2 * np.sin(theta2)
    v2 = vx2 ** 2 + vy2 ** 2

    # Compute the Potential Energy
    PE_1 = m1 * g * y1
    PE_2 = m2 * g * y2
    # PE_1 = g * m1 * (y1 + l1)
    # PE_2 = g * m2 * (y2 + l1 + l2)
    PE = PE_1 + PE_2

    # Compute the Kinetic Energy
    KE_1 = (1 / 2) * m1 * v1
    KE_2 = (1 / 2) * m2 * v2
    KE = KE_1 + KE_2

    # Compute the Total Energy
    TE = PE + KE

    return x1, y1, x2, y2, ts, theta1, theta2, omega1, omega2, PE, KE, TE


def plotAll(ts, theta1, theta2, omega1, omega2, PE, KE, TE):
    """
    Generate plots
    """
    ################################################################################
    xs = theta1
    ys = [theta2]
    labels = [""]
    xaxis = '$\\theta_1$ in $rad$'
    yaxis = '$\\theta_2$ in $rad$'
    title = '$\\theta_2$ vs $\\theta_1$'
    filename = "fig1"
    plot(xs, ys, labels, xaxis, yaxis, title, filename)
    ################################################################################
    xs = omega1
    ys = [omega2]
    labels = [""]
    xaxis = '$\omega_1$ in $rad/s$'
    yaxis = '$\omega_2$ in $rad/s$'
    title = '$\omega_2$ vs $\omega_1$'
    filename = "fig2"
    plot(xs, ys, labels, xaxis, yaxis, title, filename)
    ################################################################################
    xs = ts
    ys = [theta1, theta2]
    labels = ["$\\theta_1$", "$\\theta_2$"]
    xaxis = 'time in $s$'
    yaxis = 'Angular position in $rad$'
    title = '$\\theta_1$ and $\\theta_2$ vs time'
    filename = "fig3"
    plot(xs, ys, labels, xaxis, yaxis, title, filename)
    ################################################################################
    xs = ts
    ys = [omega1, omega2]
    labels = ["$\omega_1$", "$\omega_2$"]
    xaxis = 'time in $s$'
    yaxis = 'Angular velocity in $rad/s$'
    title = '$\omega_1$ and $\omega_2$ vs time'
    filename = "fig4"
    plot(xs, ys, labels, xaxis, yaxis, title, filename)
    ################################################################################
    xs = ts
    ys = [PE, KE, TE]
    labels = ["Potential", "Kinetic", "Total"]
    xaxis = "time"
    yaxis = "Energy in $J$"
    title = "Energies vs time"
    filename = "fig5"
    plot(xs, ys, labels, xaxis, yaxis, title, filename)
    ################################################################################
    return


def animate(x1, y1, x2, y2, ts, theta1, theta2, omega1, omega2, PE, KE, TE):
    # Misc variables
    data_skip = 60
    offset = 0.1

    # Initialize a figure
    fig = plt.figure(figsize=(15, 7))

    # Title for Initial Conditions
    units_str = r'$m_1=$' + str(m1) + 'kg, ' + \
                r'$m_2=$' + str(m2) + 'kg, ' + \
                r'$l_1=$' + str(l1) + 'm, ' + \
                r'$l_2=$' + str(l2) + 'm, ' + \
                r'$\theta_1=$' + str(theta1_0) + r'$^\circ$, ' + \
                r'$\theta_2=$' + str(theta2_0) + r'$^\circ$, ' + \
                r'$\omega_1=$' + str(omega1_0) + r'$^\circ s^{-1}$, ' + \
                r'$\omega_2=$' + str(omega2_0) + r'$^\circ s^{-1}$'

    # Display the Title of the figure
    plt.suptitle('Double Pendulum Simulation\n' +
                 'Ramsey Phuc\n' +
                 str(units_str), color='#FF00FF')

    # Initialize Double Pendulum Animation
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.grid(True)
    s = l1 + l2
    ax1.set_xlim((-s * 1.2, s * 1.2))
    ax1.set_ylim((-s * 1.2, s * 1.2))
    ax1.set_title('Double Pendulum Animation')
    plt.subplots_adjust(top=0.86, bottom=0.07)
    x_data, y_data = [], []
    pendulum1_string, = ax1.plot(x_data, y_data, lw=1, color='black')
    pendulum2_string, = ax1.plot(x_data, y_data, lw=1, color='black')
    x_path_bob2, y_path_bob2 = [], []
    trace_path_bob2, = ax1.plot(x_path_bob2, y_path_bob2, lw=1, color='black')

    # Initialize Phase Plot of bobs' position
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.set_xlabel('$\\theta_1$ in $rad$')
    ax2.set_ylabel('$\\theta_2$ in $rad$')
    ax2.set_title('$\\theta_2$ vs $\\theta_1$')
    ax2.set_xlim((min(theta1) - offset, max(theta1) + offset))
    ax2.set_ylim((min(theta2) - offset, max(theta2) + offset))
    ax2.minorticks_on()
    ax2.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    ax2.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    curve2_x, curve2_y = [], []
    curve2, = ax2.plot(curve2_x, curve2_y, lw=1, color='black')

    # Initialize Phase Plot of bobs' angular velocity
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.set_xlabel('$\omega_1$ in $rad/s$')
    ax3.set_ylabel('$\omega_2$ in $rad/s$')
    ax3.set_title('$\omega_2$ vs $\omega_1$')
    ax3.set_xlim((min(omega1) - offset, max(omega1) + offset))
    ax3.set_ylim((min(omega2) - offset, max(omega2) + offset))
    ax3.minorticks_on()
    ax3.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    ax3.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    curve3_x, curve3_y = [], []
    curve3, = ax3.plot(curve3_x, curve3_y, lw=1, color='black')

    # Initialize Time Plot of bobs' position
    ax4 = fig.add_subplot(2, 3, 4)
    ax4.set_xlabel('time in $s$')
    ax4.set_ylabel('Angular position in $rad$')
    ax4.set_title('$\\theta_1$ and $\\theta_2$ vs time')
    ax4.set_xlim((ts[0] - offset, ts[-1] + offset))
    ax4.set_ylim((min(min(theta1), min(theta2)) - offset,
                  max(max(theta1), max(theta2)) + offset))
    ax4.minorticks_on()
    ax4.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    ax4.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    curve4a_x, curve4a_y = [], []
    curve4a, = ax4.plot(curve4a_x, curve4a_y, lw=1, color='black', label=r'$\theta_1$')
    curve4b_x, curve4b_y = [], []
    curve4b, = ax4.plot(curve4b_x, curve4b_y, lw=1, color='red', label=r'$\theta_2$')
    ax4.legend(loc='best', framealpha=0.5)

    # Initialize Time Plot of bobs' angular velocity
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.set_xlabel('time in $s$')
    ax5.set_ylabel('Angular velocity in $rad/s$')
    ax5.set_title('$\omega_2$ and $\omega_1$ vs time')
    ax5.set_xlim((ts[0] - offset, ts[-1] + offset))
    ax5.set_ylim((min(min(omega1), min(omega2)) - offset,
                  max(max(omega1), max(omega2)) + offset))
    ax5.minorticks_on()
    ax5.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    ax5.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    curve5a_x, curve5a_y = [], []
    curve5a, = ax5.plot(curve5a_x, curve5a_y, lw=1, color='black', label=r'$\omega_1$')
    curve5b_x, curve5b_y = [], []
    curve5b, = ax5.plot(curve5b_x, curve5b_y, lw=1, color='red', label=r'$\omega_2$')
    ax5.legend(loc='best', framealpha=0.5)

    # Initialize Time Plot of Energies
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.set_xlabel('time in $s$')
    ax6.set_ylabel('Energy in $J$')
    ax6.set_title('Energies vs time')
    ax6.set_xlim((ts[0] - offset, ts[-1] + offset))
    ax6.set_ylim((min(min(PE), min(KE), min(TE)) - offset,
                  max(max(PE), max(KE), max(TE)) + offset))
    ax6.minorticks_on()
    ax6.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    ax6.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    curve6a_x, curve6a_y = [], []
    curve6a, = ax6.plot(curve6a_x, curve6a_y, lw=1, color='black', label='Potential')
    curve6b_x, curve6b_y = [], []
    curve6b, = ax6.plot(curve6b_x, curve6b_y, lw=1, color='red', label='Kinetic')
    curve6c_x, curve6c_y = [], []
    curve6c, = ax6.plot(curve6c_x, curve6c_y, lw=1, color='blue', label='Total')
    ax6.legend(loc='best', framealpha=0.5)

    def initialize():
        # Initialize Double Pendulum Animation
        pendulum1_string.set_data(x_data, y_data)
        pendulum2_string.set_data(x_data, y_data)
        trace_bob1 = ax1.scatter(x_data, y_data, marker='o', color='blue')
        trace_bob2 = ax1.scatter(x_data, y_data, marker='o', color='red')
        trace_path_bob2.set_data(x_path_bob2, y_path_bob2)
        bob1 = ax1.scatter((max(x1) + max(x2)), 0, marker='o', color='black')
        bob2 = ax1.scatter((max(x1) + max(x2)), -0.5, marker='o', color='red')
        bob1_text = ax1.annotate(r'$= m_1$', xy=((max(x1) + max(x2)) + 0.1, -0.04))
        bob2_text = ax1.annotate(r'$= m_2$', xy=((max(x1) + max(x2)) + 0.1, -0.55))

        # Initialize Phase Plot of bobs' position
        curve2.set_data(curve2_x, curve2_y)
        curve2_marker = ax2.scatter(curve2_x, curve2_y, color='black')

        # Initialize Phase Plot of bobs' angular velocity
        curve3.set_data(curve3_x, curve3_y)
        curve3_marker = ax3.scatter(curve3_x, curve3_y, color='black')

        # Initialize Time Plot of bobs' position
        curve4a.set_data(curve4a_x, curve4a_y)
        # curve4a_marker = ax4.scatter(curve4a_x, curve4a_y, color='black')

        curve4b.set_data(curve4b_x, curve4b_y)
        # curve4b_marker = ax4.scatter(curve4b_x, curve4b_y, color='red')

        # Initialize Time Plot of bobs' angular velocity
        curve5a.set_data(curve5a_x, curve5a_y)
        # curve5a_marker = ax5.scatter(curve5a_x, curve5a_y, color='black')

        curve5b.set_data(curve5b_x, curve5b_y)
        # curve5b_marker = ax5.scatter(curve5b_x, curve5b_y, color='red')

        # Initialize Time Plot of Energies
        curve6a.set_data(curve6a_x, curve6a_y)
        # curve6a_marker = ax6.scatter(curve6a_x, curve6a_y, color='black')

        curve6b.set_data(curve6b_x, curve6b_y)
        # curve6b_marker = ax6.scatter(curve6b_x, curve6b_y, color='red')

        curve6c.set_data(curve6c_x, curve6c_y)
        # curve6c_marker = ax6.scatter(curve6c_x, curve6c_y, color='blue')

        return [bob1, bob2, bob1_text, bob2_text,
                pendulum1_string, pendulum2_string,
                trace_bob1, trace_bob2,
                trace_path_bob2,
                curve2, curve2_marker,
                curve3, curve3_marker,
                curve4a, curve4b,
                # curve4a_marker, # curve4b_marker,
                curve5a, curve5b,
                # curve5a_marker,  # curve5b_marker,
                curve6a, curve6b, curve6c,
                # curve6a_marker,  # curve6b_marker,  # curve6c_marker,
                ]

    def update(i):
        # Updating Double Pendulum Animation
        bob1 = ax1.scatter((max(x1) + max(x2)), 0, color='black')
        bob2 = ax1.scatter((max(x1) + max(x2)), -0.5, color='red')
        bob1_text = ax1.annotate(r'$= m_1$', xy=((max(x1) + max(x2)) + 0.1, -0.04))
        bob2_text = ax1.annotate(r'$= m_2$', xy=((max(x1) + max(x2)) + 0.1, -0.55))
        pendulum1_string.set_data([0, x1[i]], [0, y1[i]])
        pendulum2_string.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        trace_bob1 = ax1.scatter(x1[i], y1[i], color='black')
        trace_bob2 = ax1.scatter(x2[i], y2[i], color='red')
        x_path_bob2.append(x2[i])
        y_path_bob2.append(y2[i])
        trace_path_bob2.set_data(x_path_bob2, y_path_bob2)

        # Updating Phase Plot of bobs' position
        curve2_x.append(theta1[i])
        curve2_y.append(theta2[i])
        curve2.set_data(curve2_x, curve2_y)
        curve2_marker = ax2.scatter(curve2_x[-1], curve2_y[-1], color='black')

        # Updating Phase Plot of bobs' angular velocity
        curve3_x.append(omega1[i])
        curve3_y.append(omega2[i])
        curve3.set_data(curve3_x, curve3_y)
        curve3_marker = ax3.scatter(curve3_x[-1], curve3_y[-1], color='black')

        # Updating Time Plot of bobs' position
        curve4a_x.append(ts[i])
        curve4a_y.append(theta1[i])
        curve4a.set_data(curve4a_x, curve4a_y)
        # curve4a_marker = ax4.scatter(curve4a_x[-1], curve4a_y[-1], color='black')

        curve4b_x.append(ts[i])
        curve4b_y.append(theta2[i])
        curve4b.set_data(curve4b_x, curve4b_y)
        # curve4b_marker = ax4.scatter(curve4b_x[-1], curve4b_y[-1], color='red')

        # Updating Time Plot of bobs' angular velocity
        curve5a_x.append(ts[i])
        curve5a_y.append(omega1[i])
        curve5a.set_data(curve5a_x, curve5a_y)
        # curve5a_marker = ax5.scatter(curve5a_x[-1], curve5a_y[-1], color='black')

        curve5b_x.append(ts[i])
        curve5b_y.append(omega2[i])
        curve5b.set_data(curve5b_x, curve5b_y)
        # curve5b_marker = ax5.scatter(curve5b_x[-1], curve5b_y[-1], color='red')

        # Updating Time Plot of Energies
        curve6a_x.append(ts[i])
        curve6a_y.append(PE[i])
        curve6a.set_data(curve6a_x, curve6a_y)
        # curve6a_marker = ax6.scatter(curve6a_x[-1], curve6a_y[-1], color='black')

        curve6b_x.append(ts[i])
        curve6b_y.append(KE[i])
        curve6b.set_data(curve6b_x, curve6b_y)
        # curve6b_marker = ax6.scatter(curve6b_x[-1], curve6b_y[-1], color='red')

        curve6c_x.append(ts[i])
        curve6c_y.append(TE[i])
        curve6c.set_data(curve6c_x, curve6c_y)
        # curve6c_marker = ax6.scatter(curve6c_x[-1], curve6c_y[-1], color='blue')

        return [bob1, bob2, bob1_text, bob2_text,
                pendulum1_string, pendulum2_string,
                trace_bob1, trace_bob2,
                trace_path_bob2,
                curve2, curve2_marker,
                curve3, curve3_marker,
                curve4a, curve4b,
                # curve4a_marker, # curve4b_marker,
                curve5a, curve5b,
                # curve5a_marker,  # curve5b_marker,
                curve6a, curve6b, curve6c,
                # curve6a_marker,  # curve6b_marker,  # curve6c_marker,
                ]

    plt.tight_layout()
    anim = FuncAnimation(fig, update, frames=np.arange(0, len(ts), data_skip),
                         init_func=initialize, interval=20, repeat=False, blit=True)
    # time1 = time.time()
    # anim.save('Double_Pendulum_Simulation.gif', dpi=300, fps=60, writer='Pillow')
    # time2 = time.time()
    # print("Save time: " + str(time2-time1) + "seconds")
    plt.show()
    return


if __name__ == '__main__':
    # Set Constants and Variables
    m1, m2 = 1, 1  # Mass [kg]
    l1, l2 = 1, 1  # Length of each string [m]
    g = 9.81  # Acceleration due to gravity [m/s^2]
    theta1_0, theta2_0 = 45, 0  # Initial Angle relative to the y-axis [deg]
    omega1_0, omega2_0 = 0, 0  # Initial angular velocity [rad/s]
    t1, t2 = 0, 10  # Time span [s]

    # Solve the system
    sol = solve(m1, m2, l1, l2, g, theta1_0, theta2_0, omega1_0, omega2_0, t1, t2)
    x1, y1, x2, y2, ts, theta1, theta2, omega1, omega2, PE, KE, TE = sol

    # Generate still plots
    plotAll(ts, theta1, theta2, omega1, omega2, PE, KE, TE)

    # Generate an animation of the double pendulum with 4 animated plots
    animate(x1, y1, x2, y2, ts, theta1, theta2, omega1, omega2, PE, KE, TE)
