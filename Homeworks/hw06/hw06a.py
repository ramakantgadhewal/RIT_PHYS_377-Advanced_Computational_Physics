"""
Author: Ramsey (Rayla) Phuc
File: hw06a.py
"""
import matplotlib.pyplot as plt
import numpy as np


def pb01():
    """
    Part a: Write a Program to compute the solution to the two-dimensional
    electrostatics problem using the Jacobi method. Take the following values:
    Box is 1 m along each sides, 𝑉 = 1 volt, and the grid spacing 𝑎 = 1 cm, so
    that there are 100 grid points on a side, or 101 if we count the points at
    both the beginning and the end.

    Part b: Make a density plot of the result.
    """
    # constants/parameters
    V = 1  # Voltage (V)

    # target accuracy
    target = 1e-6

    # Size of grid
    n = 100  # dimensions of the grid
    phi = np.zeros((n, n), dtype=float)
    phi_x = np.zeros((n, n), dtype=float)

    # Initialize potential of Grids
    for i in range(0, n):
        phi[0][i] = V
        phi_x[0][i] = V
        pass

    # Keep iterating until the target accuracy is met
    delta = 1
    while delta > target:
        for x in range(1, n - 1):
            for y in range(1, n - 1):
                top = phi[x - 1][y]
                right = phi[x][y + 1]
                bottom = phi[x + 1][y]
                left = phi[x][y - 1]
                phi_x[x][y] = (1 / 4) * (top + right + bottom + left)
                pass
            pass

        delta = np.max(abs(phi_x - phi))
        phi, phi_x = phi_x, phi

    # Create a density plot based on the data
    plt.imshow(phi)
    # Show the plot
    plt.show()
    return


def pb02():
    """
    Part a: Write a program, or modify the previous one, to solve Poisson’s equation
    of electrostatics for this system (two-dimensional), which is in presence of a
    charge density. The framework is described in the lecture.

    Part b: Make a density plot of the result.
    """
    # constants/parameters
    L = 100  # Length of the side of a square box in centimeters
    region = 20  # Size of square charge in centimeters
    epsilon_0 = 1
    a = 1 / L  # spacing between each grid point in centimeters

    # target accuracy
    target = 1e-6

    # Size of grid
    n = 100  # dimensions of the grid
    phi = np.zeros((n, n), dtype=float)
    phi_x = np.zeros((n, n), dtype=float)
    rhos = np.zeros((n, n), dtype=float)

    # Charge Density Values
    rho0, rho1, rho2 = 0, 1, -1

    # places where the potentials are initially non-zero in centimeter
    p_1 = 30, 70
    p_2 = 70, 30

    # Initialize potential of Grids
    for x in range(0, n):
        for y in range(0, n):
            # Region of Square Charge 1
            cond1r = p_1[0] - (region / 2) <= x <= p_1[0] + (region / 2)
            cond1c = p_1[1] - (region / 2) <= y <= p_1[1] + (region / 2)

            # Region of Square Charge 2
            cond2r = p_2[0] - (region / 2) <= x <= p_2[0] + (region / 2)
            cond2c = p_2[1] - (region / 2) <= y <= p_2[1] + (region / 2)
            if cond1r and cond1c:
                phi[x][y] = rho1
                phi_x[x][y] = rho1
                rhos[x][y] = rho1
            elif cond2r and cond2c:
                phi[x][y] = rho2
                phi_x[x][y] = rho2
                rhos[x][y] = rho2
            pass
        pass

    # Keep iterating until the target accuracy is met
    delta = 1
    while delta > target:
        for x in range(0, n - 1):
            for y in range(0, n - 1):
                # Boundary Condition
                cond1 = (x == 0) or (y == 0) or (x == n) or (y == n)

                # Region of Square Charge 1
                cond2r = p_1[0] - (region / 2) <= x <= p_1[0] + (region / 2)
                cond2c = p_1[1] - (region / 2) <= y <= p_1[1] + (region / 2)
                cond2 = cond2r and cond2c

                # Region of Square Charge 2
                cond3r = p_2[0] - (region / 2) <= x <= p_2[0] + (region / 2)
                cond3c = p_2[1] - (region / 2) <= y <= p_2[1] + (region / 2)
                cond3 = cond3r and cond3c

                # If the point is at the boundary or within the square charges
                if cond1 or cond2 or cond3:
                    # Keep the same value
                    phi_x[x][y] = phi[x][y]
                # Otherwise
                else:
                    # Compute new phi value
                    top = phi[x - 1][y]
                    right = phi[x][y + 1]
                    bottom = phi[x + 1][y]
                    left = phi[x][y - 1]
                    term = rhos[x][y] * a ** 2 / (4 * epsilon_0)
                    phi_x[x][y] = (1 / 4) * (top + right + bottom + left) + term
                pass
            pass

        delta = np.max(abs(phi_x - phi))
        phi, phi_x = phi_x, phi

    # Create a density plot based on the data
    plt.imshow(phi, cmap="hot")
    # Show the plot
    plt.show()
    return


if __name__ == '__main__':
    pb01()
    pb02()
