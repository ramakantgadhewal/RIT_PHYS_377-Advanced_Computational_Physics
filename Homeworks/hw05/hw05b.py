"""
Author: Ramsey (Rayla) Phuc
File: hw05b.py
"""
import matplotlib.pyplot as plt
import numpy as np


##########################################################################################
# Helper Functions
##########################################################################################
def get_data(data):
    """
    "Flattens" the list. Converts a (n x m) array to a (m x n) list
    """
    numRows, numCols = len(data), len(data[0])
    temp = np.zeros((numCols, numRows), dtype=float)
    for r in range(0, numRows):
        for c in range(0, numCols):
            temp[c][r] = data[r][c]
    return np.asarray(temp, dtype=float)


def rk4(df, r, ts, h, E):
    """
    Runge-Kutta 4 Method

    :param E:
    :param df: System of ODE(s) to solve
    :param r: Initial Condition vector
    :param ts: time
    :param h: Step Size
    """
    # Initialize the approximated solution
    rs = [r]

    # For each x value
    for i in range(1, len(ts)):
        # Compute f1 = f(r(t), t)
        f1 = df(rs[i - 1], ts[i - 1], E)

        # Compute f2 = f(r(t) + (h/2)*f1, t + (h/2))
        f2 = df(rs[i - 1] + (h / 2) * f1, ts[i - 1] + (h / 2), E)

        # Compute f3 = f(r(t) + (h/2)*f2, t + (h/2))
        f3 = df(rs[i - 1] + (h / 2) * f2, ts[i - 1] + (h / 2), E)

        # Compute f4 = f(r(t) + h*f3, t + h)
        f4 = df(rs[i - 1] + h * f3, ts[i - 1] + h, E)

        # Compute r(t + h) = r(t) + (h/6)*(f1 + 2*f2 + 2*f3 + f4)
        r = rs[i - 1] + (h / 6) * (f1 + 2 * f2 + 2 * f3 + f4)

        # Record the approximated value
        rs.append(r)

    # Return rs as an array
    return get_data(np.asarray(rs, dtype=float))


##########################################################################################
# Main Functions
##########################################################################################
def solve(df, r, ts, h, e, Es):
    E1, E2 = Es
    # Solve the Equation and get the right bound value
    psi2 = rk4(df, r, ts, h, E1)[0][-1]

    # Do Secant method until you get the correct energy
    target = e / 1000
    while abs(E1 - E2) > target:
        sol = rk4(df, r, ts, h, E2)[0]
        psi1, psi2 = psi2, sol[-1]
        E1, E2 = E2, E2 - psi2 * (E2 - E1) / (psi2 - psi1)
    return sol, E2 / e


def pb01():
    def V1(x):
        """
        Part a
        Potential Energy at position x
        """
        return V0 * (x / a) ** 2

    def df1(F, x, E):
        """
        Part a
        Time-Independent Schrodinger Equation in 1-D
        """
        psi, phi = F
        psi_x = phi
        phi_x = (2 * m / hbar ** 2) * (V1(x) - E) * psi
        return np.array([psi_x, phi_x], dtype=float)

    def V2(x):
        """
        Part b
        Potential Energy at position x
        """
        return V0 * (x / a) ** 4

    def df2(F, x, E):
        """
        Part b
        Time-Independent Schrodinger Equation in 1-D
        """
        psi, phi = F
        psi_x = phi
        phi_x = (2 * m / hbar ** 2) * (V2(x) - E) * psi
        return np.array([psi_x, phi_x], dtype=float)

    # Constants
    hbar = 1.0546e-34  # kg ⋅ m^2 ⋅ s^(-1)
    e = 1.6022e-19  # e

    # Parameters
    m = 9.1094e-31  # kg
    # V0 = 0 * e  # eV
    V0 = 50 * e  # eV
    # a = 5.2918e-11  # m
    a = 1e-11  # m
    w = 10 * a
    E1, E2 = 100 * e, 200 * e  # eV

    # Initial Conditions
    psi0, phi0 = 0, 1
    IC = np.array([psi0, phi0], dtype=float)

    # Time Interval
    t1, t2 = -w, w

    # Number of intervals
    n = 1000

    # Step Size
    h = 2 * w / n

    # t values
    ts = np.linspace(t1, t2, n + 1)

    print("Part a: ")
    # Solve for Ground State Energy
    sol = solve(df1, IC, ts, h, e, [E1, E2])
    print("The calculated energy at the ground state is " + str(sol[1]) + " eV")

    # Solve for 1st Excited State Energy
    E1, E2 = 300 * e, 310 * e  # eV
    sol = solve(df1, IC, ts, h, e, [E1, E2])
    print("The calculated energy at the first excited state is " + str(sol[1]) + " eV")

    # Solve for 2nd Excited State Energy
    E1, E2 = 600 * e, 610 * e
    sol = solve(df1, IC, ts, h, e, [E1, E2])
    print("The calculated energy at the second excited state is " + str(sol[1]) + " eV")

    print("\nPart b: ")
    # Solve for Ground State Energy
    E1, E2 = 100 * e, 200 * e  # eV
    sol = solve(df2, IC, ts, h, e, [E1, E2])
    print("The calculated energy at the ground state is " + str(sol[1]) + " eV")

    # Solve for 1st Excited State Energy
    E1, E2 = 300 * e, 310 * e  # eV
    sol = solve(df2, IC, ts, h, e, [E1, E2])
    print("The calculated energy at the first excited state is " + str(sol[1]) + " eV")

    # Solve for 2nd Excited State Energy
    E1, E2 = 600 * e, 610 * e
    sol = solve(df2, IC, ts, h, e, [E1, E2])
    print("The calculated energy at the second excited state is " + str(sol[1]) + " eV")

    print("\nPart c:")
    # Calculating the normalized wavefunctions of the anharmonic oscilator for the
    # ground, 1st excited, and 2nd excited states.

    # Changing width from 10a t0 5a
    w = 5 * a
    t1, t2 = -w, w
    h = 2 * w / n
    ts = np.linspace(t1, t2, n + 1)

    # Store solutions and labels for part c
    sols, labels = [], ["Ground State", "1st Excited State", "2nd Excited State"]

    E1, E2 = 100 * e, 200 * e  # eV
    sol = solve(df2, IC, ts, h, e, [E1, E2])
    # print("The calculated energy at the ground state is " + str(sol[1]) + " eV")
    sols.append(sol[0])

    # Solve for 1st Excited State Energy
    E1, E2 = 300 * e, 310 * e  # eV
    sol = solve(df2, IC, ts, h, e, [E1, E2])
    # print("The calculated energy at the first excited state is " + str(sol[1]) + " eV")
    sols.append(sol[0])

    # Solve for 2nd Excited State Energy
    E1, E2 = 600 * e, 610 * e
    sol = solve(df2, IC, ts, h, e, [E1, E2])
    # print("The calculated energy at the second excited state is " + str(sol[1]) + " eV")
    sols.append(sol[0])

    # For each solution
    for i in range(0, len(sols)):
        sol = sols[i]
        # Compute the area under the curve
        res = 0
        for j in range(0, len(sol) - 1):
            res += 0.5 * (np.square(sol[j]) + np.square(sol[j + 1]))
        # Plot the normalized wave function
        sol /= np.sqrt(res)
        plt.plot(ts, sol, label=labels[i])
        pass

    plt.legend()
    plt.show()
    return 0


if __name__ == '__main__':
    pb01()
