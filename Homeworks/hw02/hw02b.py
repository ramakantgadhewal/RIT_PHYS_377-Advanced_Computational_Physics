"""
Author: Ramsey (Rayla) Phuc
File: hw02b.py
"""

# Imported to use complex mathematical operations
import math

# Imported to create and display a plot
import matplotlib.pyplot as plt
# Imported to create evenly spaced numbers
import numpy as np


def simpson(f, a, b, n):
    """
    Performs the simpson's rule for approximating area
    :param f: f(x)
    :param a: lower bound
    :param b: upper bound
    :param n: Number of trapezoids
    :return: The area under f(x)
    """
    # Get the interval step
    h = (b - a) / n

    # Add the value of f at the endpoints
    res = f(a) + f(b)

    # Doing summations for the odd case
    for i in range(1, n - 1, 2):
        res += 4 * f(a + i * h)
        pass

    # Doing summations for the even case
    for j in range(2, n - 2, 2):
        res += 2 * f(a + j * h)
        pass

    # multiply the current sum by h/3
    res *= h / 3

    # Return the result
    return res


def pb01():
    print("Part a: Plotting Bessel functions.")
    a, b, n = 0, math.pi, 1000

    def J(m, x):
        """
        Calculates the value of the Bessel function using Simpson's rule with
        N = 1000 points.
        """

        def f(theta):
            """
            Integrand of the Bessel Function
            """
            return math.cos(m * theta - x * math.sin(theta))

        return (1 / math.pi) * simpson(f, a, b, n)

    # Setting up x values
    x1, x2 = 0, 20

    # Initialize x values
    xs = np.linspace(x1, x2, n)

    # For each Bessel function of the following m values
    ms = [0, 1, 2]
    for m in ms:
        # Get the y values at each x value
        ys = []
        for x in xs:
            ys.append(J(m, x))
            pass
        # Plot x vs y
        plt.plot(xs, ys, label="J_" + str(m) + "(x)")
        pass

    # Display the legend
    plt.legend()

    # Display the gridlines
    plt.minorticks_on()
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='black')
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

    # Display the plot
    plt.show()

    print("\nPart b: Creating a density plot of the intensity of the circular "
          "\ndiffraction pattern.")

    def I(r):
        """
        The equation that represents the intensity of the light in the
        diffraction pattern
        """
        return ((J(1, k * r)) / (k * r)) ** 2

    # Point light source in nanometers
    wavelength = 500

    # Converting lambda to k
    k = (2 * math.pi) / wavelength

    factor = 2

    # Distance in the focal plane from the center of the diffraction pattern
    # in nanometers
    r1, r2 = 0, 1 * factor * (10 ** 3)

    # Center of the plot
    c = r2 / 2

    # Step size
    q = 5 * factor

    # Creating data for density plot
    data = []
    for i in range(r1, r2, q):
        print(i)
        temp = []
        for j in range(r1, r2, q):
            # Calculating the distance from the center to a point
            r = math.sqrt((i - c) ** 2 + (j - c) ** 2)

            # Calculate the intensity given the distance
            if k * r == 0:
                # Special case where k*r = 0
                intensity = (1 / 2) ** 2
            else:
                intensity = I(r)

            # Add the calculated intensity to the data
            temp.append(intensity)
            pass

        data.append(temp)
        pass

    # Create a density plot based on the data
    plt.imshow(data, vmax=0.01, cmap="hot")

    # Show the plot
    plt.show()
    pass


if __name__ == '__main__':
    pb01()
    pass
