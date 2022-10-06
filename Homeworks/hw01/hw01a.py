"""
Author: Ramsey (Rayla) Phuc
File: hw01a.py
"""

# Imported to use complex mathematical operations
import math

# Imported to create and display a plot
import matplotlib.pyplot as plt
# Imported to read a .txt file
import numpy as np


def pb01():
    """
    A ball is dropped from a tower of height height with initial velocity zero.
    Write a program that asks the user to enter the height of the tower in meters
    and then calculates and prints the time the ball takes until it hits the
    ground, ignoring air resistance. Use your program to calculate the time for
    a ball dropped from a 100 m high tower.
    """
    # Gets the user inputted height
    height = float(input("Enter the height of the tower in meters: "))
    # Computes the time it takes for the object to hit the ground
    # The formula used is one of the four Kinematics equations
    t = math.sqrt(height * 2 / 9.81)
    # Prints the results
    print("For height=" + str(height) + ", t=" + str(t))

    # Set the height to 100
    height = 100
    # Computes the time it takes for the object to hit the ground
    # The formula used is one of the four Kinematics equations
    t = math.sqrt(height * 2 / 9.81)
    # Prints the results
    print("For height=" + str(height) + ", t=" + str(t))
    pass


def pb02():
    """
    Suppose the position of a point in two-dimensional space is given to us in
    polar coordinates r, θ. We want to convert it to Cartesian coordinates x, y.
    Write a program to do this? The appropriate steps are: (1) Get the user to
    enter the values of r and θ (2) Convert those values to cartesian coordinates
    using the standard formulas: x = rcosθ, y = rsinθ (3) Print out the results.
    """
    # Gets the Polar coordinates from user input
    r = float(input("Enter r: "))
    theta = float(input("Enter θ in radians: "))
    # Convert from Polar to Cartesian coordinates
    x, y = r * math.cos(theta), r * math.sin(theta)
    # Print the Cartesian coordinates
    print("x=" + str(x) + ", y=" + str(y))
    pass


def pb03():
    """
    Suppose the position of a point in two-dimensional space is given to us in
    Cartesian coordinates x, y. We want to convert it to polar coordinates r, θ.
    Write a program to do this? The appropriate steps are: (1) Get the user to
    enter the values of x and y (2) Convert those values to polar coordinates
    using the standard formulas: x = rcosθ, y = rsinθ (3) Print out the results.
    """
    # Gets the Cartesian coordinates from user input
    x = float(input("Enter x: "))
    y = float(input("Enter y: "))
    # Convert from Cartesian to Polar coordinates
    if (x > 0) and (y > 0):
        theta = math.atan(y / x)
    elif (x < 0) and (y > 0):
        theta = 90 + math.atan(y / x)
    elif (x < 0) and (y < 0):
        theta = 180 + math.atan(y / x)
    else:
        theta = 360 - math.atan(y / x)
    r = x / math.cos(theta)
    # Print the Polar coordinates
    print("r=" + str(r) + ", θ=" + str(theta))
    pass


def pb04():
    """
    Calculate the distance of a point in cylindrical coordinates to the origin.
    Print out the results.
    """
    # Gets the Cylindrical coordinates from user input
    rho = float(input("Enter ρ: "))
    theta = float(input("Enter θ in radians: "))
    z = float(input("Enter z: "))
    # Computes the distance between the point and the origin in Cylindrical
    # coordinates
    d = math.sqrt(rho ** 2 - 2 * rho * math.cos(theta) + z ** 2)
    # Print the results
    print("d=" + str(d))
    pass


def pb05():
    """
    There is a file on myCourses folder “Files and Data” called stm.txt, which
    contains a grid of values from scanning tunneling microscope measurements
    of the (111) surface of silicon. A scanning tunneling microscope (STM) is a
    device that measures the shape of a surface at the atomic level by tracking
    a sharp tip over the surface and measuring quantum tunneling current as a
    function of position. The end result is a grid of values that represent the
    height of the surface and the file stm.txt contains just such a grid of
    values.

    Write a program that reads the data contained in the file and makes a
    density plot of the values. Use the various options and variants you have
    learned about to make a picture that shows the structure of the silicon
    surface clearly.

    *** Need to have stm.txt in the same location as hw01a.py***
    """
    # Reads data
    data = np.loadtxt('stm.txt')
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
    print("\n-----Problem 4-----")
    # pb04()
    print("\n-----Problem 5-----")
    # pb05()
    pass
