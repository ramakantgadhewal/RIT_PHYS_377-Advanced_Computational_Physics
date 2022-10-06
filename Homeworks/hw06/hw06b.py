"""
Author: Ramsey (Rayla) Phuc
File: hw06b.py
"""
import matplotlib.pyplot as plt
import numpy as np


def pb01(omega):
    """
    Write a Program that implements the combined overrelaxation/Gauss–Seidel method to
    solve Laplace’s equation for the two-dimensional problem
    (similar to HW Problem 6a-1)—a square box 1 m on each side, at voltage V = 1 volt
    along the top wall and zero volts along the other three. Use a grid of spacing
    a = 1 cm, so that there are 100 grid points along each wall, or 101 if you count
    the points at both ends. Continue the iteration of the method until the value of
    the electric potential changes by no more than δ = 10−6 V at any grid point on any
    step, then make a density plot of the final solution. Experiment with different
    values of ω to find which value gives the fastest solution. In general, larger
    values cause the calculation to run faster, but if you choose too large a value
    the speed drops off and for values above 1 the calculation becomes unstable.
    """
    # constants/parameters
    V = 1  # Voltage (V)

    # target accuracy
    target = 1e-6

    # Size of grid
    n = 100  # dimensions of the grid
    phi = np.zeros((n, n), dtype=float)

    # Initialize potential of Grids
    for i in range(0, n):
        phi[0][i] = V
        pass

    # Keep iterating until the target accuracy is met
    delta = 1
    print(delta)
    while delta > target:
        delta = 0
        for x in range(1, n - 1):
            for y in range(1, n - 1):
                top = phi[x - 1][y]
                right = phi[x][y + 1]
                bottom = phi[x + 1][y]
                left = phi[x][y - 1]
                term = omega * phi[x][y]
                phi_new = ((1 + omega) / 4) * (top + right + bottom + left) - term
                delta_n = abs(phi_new - phi[x][y])
                phi[x][y] = phi_new
                if delta < delta_n:
                    delta = delta_n
                    pass
                pass
            pass
        print(delta)
        pass

    # Create a density plot based on the data
    plt.imshow(phi)
    # Show the plot
    plt.show()
    return


def pb02(omega):
    """
    Using any of the methods we have studied, write a program to calculate the
    electrostatic potential in the box on a grid of 100 × 100 points, where the
    walls of the box are at voltage zero and the two plates (which are of negligible
    thickness) are at voltages ±1 V as shown. Have your program calculate the value
    of the potential at each grid point to a precision of 10−6 volts and then make a
    density plot of the result.
    """
    # constants/parameters
    L = 100  # Length of the side of a square box in centimeters
    plate_l = 60  # size of plate in centimeters

    # target accuracy
    target = 1e-6

    # Size of grid
    n = 100  # dimensions of the grid
    phi = np.zeros((n, n), dtype=float)
    phi_x = np.zeros((n, n), dtype=float)

    # Initial Voltages
    V1, V2 = 1, -1

    # places where the plates are initially non-zero in centimeters
    plate1 = 50, 20
    plate2 = 50, 80

    # Initialize potential of Grids
    for r in range(0, n):
        for c in range(0, n):
            # Region of plate 1
            cond1r = plate1[0] - (plate_l / 2) <= r <= plate1[0] + (plate_l / 2)
            cond1c = c == plate1[1]

            # Region of plate 2
            cond2r = plate2[0] - (plate_l / 2) <= r <= plate2[0] + (plate_l / 2)
            cond2c = c == plate2[1]

            if cond1r and cond1c:
                phi[r][c] = V1
                phi_x[r][c] = V1
                pass
            if cond2r and cond2c:
                phi[r][c] = V2
                phi_x[r][c] = V2
                pass
            pass
        pass

    # Keep iterating until the target accuracy is met
    delta = 1
    while delta > target:
        # delta = 0
        print(delta)
        for x in range(0, n - 1):
            for y in range(0, n - 1):
                # Boundary Condition
                cond = (x == 0) or (y == 0) or (x == n) or (y == n)

                # Region of plate 1
                cond1r = plate1[0] - (plate_l / 2) <= x <= plate1[0] + (plate_l / 2)
                cond1c = y == plate1[1]
                cond1 = cond1r and cond1c

                # Region of plate 2
                cond2r = plate2[0] - (plate_l / 2) <= x <= plate2[0] + (plate_l / 2)
                cond2c = y == plate2[1]
                cond2 = cond2r and cond2c

                # If the point is at the boundary or on the plate
                if cond1 or cond2 or cond:
                    # Keep the same value
                    phi_x[x][y] = phi[x][y]
                    pass
                # Otherwise
                else:
                    # Compute new phi value
                    top = phi[x - 1][y]
                    right = phi[x][y + 1]
                    bottom = phi[x + 1][y]
                    left = phi[x][y - 1]
                    phi_x[x][y] = (1 / 4) * (top + right + bottom + left)
                    # phi_new = (1 / 4) * (top + right + bottom + left)
                    # delta_n = abs(phi_new - phi[x][y])
                    # phi[x][y] = phi_new
                    # if delta < delta_n:
                    #     delta = delta_n
                    #     pass
                    pass
                pass
            pass
        delta = np.max(abs(phi_x - phi))
        phi, phi_x = phi_x, phi
        pass

    # Create a density plot based on the data
    plt.imshow(phi, cmap="hot")
    # Show the plot
    plt.show()
    return


def pb02_test(omega):
    """
    Using any of the methods we have studied, write a program to calculate the
    electrostatic potential in the box on a grid of 100 × 100 points, where the
    walls of the box are at voltage zero and the two plates (which are of negligible
    thickness) are at voltages ±1 V as shown. Have your program calculate the value
    of the potential at each grid point to a precision of 10−6 volts and then make a
    density plot of the result.
    """
    # constants/parameters
    L = 100  # Length of the side of a square box in centimeters
    plate_l = 60  # size of plate in centimeters

    # Initial Voltages
    V1, V2 = 1, -1

    # target accuracy
    target = 1e-6

    # Size of grid
    n = 100  # dimensions of the grid
    phi = np.zeros((n, n), dtype=float)

    # places where the plates are initially non-zero in centimeters
    plate1 = 50, 20
    plate2 = 50, 80

    # Initialize potential of Grids
    for r in range(0, n):
        for c in range(0, n):
            # Region of plate 1
            cond1r = plate1[0] - (plate_l / 2) <= r <= plate1[0] + (plate_l / 2)
            cond1c = c == plate1[1]

            # Region of plate 2
            cond2r = plate2[0] - (plate_l / 2) <= r <= plate2[0] + (plate_l / 2)
            cond2c = c == plate2[1]

            if cond1r and cond1c:
                phi[r][c] = V1
                pass
            if cond2r and cond2c:
                phi[r][c] = V2
                pass
            pass
        pass

    # Keep iterating until the target accuracy is met
    delta = 1
    while delta > target:
        delta = 0
        for x in range(0, n - 1):
            for y in range(0, n - 1):
                # Boundary Condition
                cond = (x == 0) or (y == 0) or (x == n) or (y == n)

                # Region of plate 1
                cond1r = plate1[0] - (plate_l / 2) <= x <= plate1[0] + (plate_l / 2)
                cond1c = y == plate1[1]
                cond1 = cond1r and cond1c

                # Region of plate 2
                cond2r = plate2[0] - (plate_l / 2) <= x <= plate2[0] + (plate_l / 2)
                cond2c = y == plate2[1]
                cond2 = cond2r and cond2c

                # If the point is not at the boundary or on the plates
                if not (cond1 or cond2 or cond):
                    # Compute new phi value
                    top = phi[x - 1][y]
                    right = phi[x][y + 1]
                    bottom = phi[x + 1][y]
                    left = phi[x][y - 1]
                    term = omega * phi[x][y]
                    phi_new = ((1 + omega) / 4) * (top + right + bottom + left) - term
                    delta_n = abs(phi_new - phi[x][y])
                    phi[x][y] = phi_new
                    if delta < delta_n:
                        delta = delta_n
                        pass
                    pass
                pass
            pass
        pass

    # Create a density plot based on the data
    plt.imshow(phi, cmap="hot")
    # Show the plot
    plt.show()
    return


if __name__ == '__main__':
    # times = []
    # omegas = np.linspace(0.936, 0.938, 3)
    # for omega in omegas:
    #     start = time.time()
    #     pb01(omega)
    #     end = time.time()
    #     times.append(abs(start - end))
    # for i in range(0, len(times)):
    #     print("Omega = " + str(omegas[i]) + ": " + str(times[i]) + " seconds")
    #     pass

    """
    On my laptop, the value of omega that performed the algorithm the fastest was omega = 0.938
    """
    omega = 0.938
    pb01(omega)
    # omega = 0.5
    pb02(omega)
    # pb02_test(omega)
