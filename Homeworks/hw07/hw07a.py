"""
Author: Ramsey (Rayla) Phuc
File: hw07a.py
"""
import matplotlib.pyplot as plt
import numpy as np


def pb01():
    """
    Write a program, or modify an earlier program, to calculate the temperature
    profile of the crust as a function of depth up to 20 m and time up to 10 years.
    Start with temperature everywhere equal to 10°C, except at the surface and the
    deepest point, choose values for the number of grid points and the time-step h,
    then run your program for the first nine simulated years, to allow it to settle
    down into whatever pattern it reaches. Then for the tenth and final year plot
    four temperature profiles taken at 3-month intervals on a single graph to
    illustrate how the temperature changes as a function of depth and time.
    """

    def T(t):
        """
        Function to represent the mean daily temperature at a particular point on
        the surface
        """
        return A + B * np.sin(2 * np.pi * t / tau)

    # Constants and Parameters
    A = 10  # Temperature [C]
    B = 12  # Temperature [C]
    tau = 365  # Time in [days]
    L = 20  # [m]
    D = 0.1  # Thermal diffusivity [m^2 day^(-1)]

    # Different temperatures
    TLow = 11
    TMid = 10
    THi = 12

    # Time interval
    tMin, tMax = 9 * 365, 10 * 365  # Time interval in [days]
    ts = np.linspace(tMin, tMax, 4 + 1)

    # Initialize 2 arrays with size n
    n = 100
    T1, T2 = np.zeros(n, dtype=float), np.zeros(n, dtype=float)

    # Initializing the arrays
    T1[0], T1[1:n - 1], T1[n - 1] = THi, TMid, TLow
    T2[0], T2[n - 1] = THi, TLow

    # Grid Spacing
    a = L / n

    # Time step
    h = 0.01
    epsilon = h / 1000

    # Conditions to end the loop
    t = 0
    tEnd = tMax + epsilon
    while t < tEnd:
        print(t / 365)
        # Compute new temperature values
        T2[0] = T(t)
        for i in range(1, n - 1):
            T2[i] = T1[i] + (h * D / (a ** 2)) * (T1[i + 1] + T1[i - 1] - 2 * T1[i])
        # Swap arrays
        T1, T2 = T2, T1
        # Increment time counter
        t += h
        # Plot the 10th year in 3-month intervals
        for tVal in ts[1:]:
            if abs(t - tVal) < epsilon:
                plt.plot(T1, label=str(tVal / 365))
    plt.legend()
    plt.savefig("hw07a.png")
    plt.show()

    return


if __name__ == '__main__':
    pb01()
