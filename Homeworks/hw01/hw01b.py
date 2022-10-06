"""
Author: Ramsey (Rayla) Phuc
File: hw01b.py
"""

# Imported to use complex mathematical operations
import math

# Imported to create and display a plot
import matplotlib.pyplot as plt


def pb01():
    """
    Quadratic Equations problem
    """

    def quadratic_v1(a, b, c):
        """
        Computes the roots of a quadratic equation given the coefficients a, b,
        and c using the standard formula.
        :param a: Coefficient of x^2
        :param b: Coefficient of x
        :param c: Constant Coefficient
        :return: Roots of a quadratic equation given the coefficients a, b, and
                 c using the standard formula.
        """
        x1 = (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
        x2 = (-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
        return x1, x2

    def quadratic_v2(a, b, c):
        """
        Computes the roots of a quadratic equation given the coefficients a, b,
        and c using the alternative formula.
        :param a: Coefficient of x^2
        :param b: Coefficient of x
        :param c: Constant Coefficient
        :return: Roots of a quadratic equation given the coefficients a, b, and
                 c using the alternative formula.
        """
        x1 = (2 * c) / (-b - math.sqrt(b ** 2 - 4 * a * c))
        x2 = (2 * c) / (-b + math.sqrt(b ** 2 - 4 * a * c))
        return x1, x2

    # Computing the roots of the quadratic equation with
    # a = 0.001, b = 1000, c = 0.001
    roots = quadratic_v1(0.001, 1000, 0.001)
    # Printing the results
    print("x1=" + str(roots[0]) + ", x2=" + str(roots[1]))

    # Computing the roots of the quadratic equation with
    # a = 0.001, b = 1000, c = 0.001
    roots = quadratic_v2(0.001, 1000, 0.001)
    # Printing the results
    print("x1=" + str(roots[0]) + ", x2=" + str(roots[1]))

    """
    Both formulas of the quadratic equation yield slightly similar results. The
    roots are not equivalent, but rather they are off by a small margin.  
    """
    pass


def pb02():
    """
    Calculating Derivatives problem
    """

    def f(x):
        """
        Function to evaluate
        :param x: variable x
        :return: f(x)
        """
        return x * (x - 1)

    def derivative(f, x, delta):
        """
        Computes the derivative of f at a point x with a small change in delta
        :param f: Function to take the derivative of
        :param x: variable x
        :param delta: small change
        :return: Derivative of f at a point x with a small change in delta
        """
        return (f(x + delta) - f(x)) / delta
        pass

    # Computing the derivative of f given x and delta
    x, delta = 1, 10 ** -2
    res = derivative(f, x, delta)
    # Print the results
    print("delta = " + str(delta) + ", df/dx = " + str(res))

    """
    The result does not perfectly match with the analytical answer, which is 1. 
    The reason why is because our delta is not small enough.
    """

    # Different deltas to test
    deltas = [10 ** -4, 10 ** -6, 10 ** -8, 10 ** -10, 10 ** -12, 10 ** -14]
    # Test each delta with a for loop
    for delta in deltas:
        # Computing the derivative of f given x and delta
        res = derivative(f, x, delta)
        # Print the results
        print("delta = " + str(delta) + ", df/dx = " + str(res))

    """
    The reason why the derivative gets more inaccurate as delta gets even closer 
    to 0 is due to arithmetic underflow. 
    """
    pass


def pb03():
    """
    Wave Interference problem
    """

    def radius(x, y, xc, yc):
        """
        Calculates the distance from a point to the center of a circle
        :param x: x coordinate of a point
        :param y: y coordinate of a point
        :param xc: x coordinate of circle center
        :param yc: y coordinate of circle center
        :return:
        """
        return math.sqrt((x - xc) ** 2 + (y - yc) ** 2)

    def height(a0, k, r1, r2):
        """
        Calculates the height of a sine wave
        :param a0: Amplitude of the waves
        :param k: Wavevector
        :param r1: Distance from a point to the center of circle 1
        :param r2: Distance from a point to the center of circle 2
        :return: Height of a sine wave
        """
        return a0 * (math.sin(k * r1) + math.sin(k * r2))

    def getK(wavelength):
        """
        Wavelength to wavevector conversion
        :param wavelength: value of wavelength
        :return: value of wavevector
        """
        return 2 * math.pi / wavelength

    # Initial variables
    wavelength, amplitude = 5, 1  # in cm

    # Calculating centers
    x1, y1 = 200 * 0.2, 250 * 0.2
    x2, y2 = 300 * 0.2, 250 * 0.2

    # Creating data for density plot
    data = []
    for i in range(0, 500):
        temp = []
        for j in range(0, 500):
            # Compute the radius of each circle at (i, j)
            r1 = radius(i * 0.2, j * 0.2, x1, y1)
            r2 = radius(i * 0.2, j * 0.2, x2, y2)
            # Compute the height of the wave
            h = height(amplitude, getK(wavelength), r1, r2)
            # Add the value of the height of the wave to the data
            temp.append(h)
        data.append(temp)

    # Create a density plot based on the data
    plt.imshow(data)

    # Show the plot
    plt.show()
    pass


if __name__ == '__main__':
    """
    Uncomment the line to run the code for that particular problem. 
    """
    print("-----Problem 1-----")
    # pb01()
    print("\n-----Problem 2-----")
    # pb02()
    print("\n-----Problem 3-----")
    pb03()
    pass
