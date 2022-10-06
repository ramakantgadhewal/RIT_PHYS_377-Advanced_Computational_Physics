"""
Author: Ramsey (Rayla) Phuc
"""
import random

import matplotlib.pyplot as plt
import numpy as np
from vpython import *


def pb01():
    """
    Code for Problem 01
    """

    """
    Part a: Write a program that generates and prints out two random numbers between
    1 and 6, to simulate the rolling of two dice.
    """
    num1, num2 = random.randrange(1, 7, 1), random.randrange(1, 7, 1)
    print("Die 1: " + str(num1) + ", Die 2: " + str(num2))

    """
    Part b: Modify your program to simulate the rolling of two dice a million times 
    and count the number of times you get a double six. Divide by a million to get 
    the fraction of times you get a double six. You should get something close to, 
    though probably not exactly equal to 1/36.
    """
    count = 0
    n = 1000000
    for _ in range(0, n):
        num1, num2 = random.randrange(1, 7, 1), random.randrange(1, 7, 1)
        if num1 == num2 == 6:
            count += 1
    print("Chances of Double 6's: " + str(count) + " / " + str(n) + " = " + str(count / n))
    return


def pb02():
    """
    Code for Problem 02
    """

    def P(t, tau):
        """
        The probability for any particular atom to decay in a time interval length t
        """
        return 1 - 2 ** (-t / tau)

    """
    Part a: For each atom of 209Pb in turn, decide at random, with the appropriate 
    probability, whether it decays or not. The probability can be calculated from 
    (class notes) is P(t). Count the total number that decay, subtract it from the 
    number of 209Pb atoms, and add it to the number of 209Bi atoms. 
    """

    # Constants and Parameters
    Bi_213 = 10000  # Number of Bi-213 atoms
    Ti_209 = 0  # Number of Ti-209 atoms
    Pb_209 = 0  # Number of Pb-209 atoms
    Bi_209 = 0  # Number of Bi-209 atoms
    tau_Bi_213 = 46 * 60  # Half-life of Bi-213
    tau_Ti_209 = 2.2 * 60  # Half-life of Ti-209
    tau_Pb_209 = 3.3 * 60  # Half-life of Pb-209

    # Time interval
    t1, t2 = 0, 20000
    h = 1  # Time step
    ts = np.arange(t1, t2, h)

    # Initialize the data set for each atom
    Bi_213_points = []
    Ti_209_points = []
    Pb_209_points = []
    Bi_209_points = []

    # For each time step
    for _ in ts:
        # Record the number of atoms for each type of atom
        Bi_213_points.append(Bi_213)
        Ti_209_points.append(Ti_209)
        Pb_209_points.append(Pb_209)
        Bi_209_points.append(Bi_209)

        decay = 0

        """
        Part a: For each atom of 209Pb in turn, decide at random, with the appropriate 
        probability, whether it decays or not. The probability can be calculated from 
        P(t, tau). Count the total number that decay, subtract it from the number of 
        209Pb atoms, and add it to the number of 209Bi atoms.
        """
        for _ in range(0, Pb_209):
            if random.random() < P(h, tau_Pb_209):
                decay += 1
        Pb_209 -= decay
        Bi_209 += decay

        """
        Part b: Now do the same for 209Tl, except that decaying atoms are subtracted 
        from the total for 209Tl and added to the total for 209Pb.
        """
        for _ in range(0, Ti_209):
            if random.random() < P(h, tau_Ti_209):
                decay += 1
        Ti_209 -= decay
        Pb_209 += decay

        """
        Part c: For 213Bi the situation is more complicated: when a 213Bi atom decays 
        you have to decide at random with the appropriate probability the route by which 
        it decays. Count the numbers that decay by each route and add and subtract 
        accordingly.
        """
        for _ in range(0, Bi_213):
            if random.random() < P(h, tau_Bi_213):
                decay += 1
        Bi_213 -= decay
        rand = random.random() * 100
        if 0 <= rand <= 97.91:
            Pb_209 += decay
        else:
            Ti_209 += decay

    plt.plot(ts, Bi_213_points, label="Bi-213")
    plt.plot(ts, Ti_209_points, label="Ti-209")
    plt.plot(ts, Pb_209_points, label="Pb-209")
    plt.plot(ts, Bi_209_points, label="Bi-209")
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Number of atoms")
    plt.savefig("hw08a_pb01.png")
    plt.show()
    return


def pb03():
    """
    Code for Problem 03
    """

    """
    Write a program to perform a million steps of this process on a lattice with 
    L = 101 and make an animation on the screen of the position of the particle. 
    (We choose an odd length for the side of the square so that there is one lattice 
    site exactly in the center.)
    """

    # Set the dimensions of the box
    side = 101
    thk = 0.1
    s2 = 2 * side - thk
    s3 = 2 * side + thk
    n = 1000000

    mid = (side - 1) / 2

    # Create the walls
    wallR = box(pos=vector(side, 0, 0), size=vector(thk, s2, s3))
    wallL = box(pos=vector(-side, 0, 0), size=vector(thk, s2, s3))
    wallB = box(pos=vector(0, -side, 0), size=vector(s3, thk, s3))
    wallT = box(pos=vector(0, side, 0), size=vector(s3, thk, s3))
    wallBK = box(pos=vector(0, 0, -side), size=vector(s2, s2, thk))

    # Create a sphere to represent a particle
    ball = sphere(color=color.green, radius=0.5, make_trail=True, retain=200)
    # ball.mass = 1.0
    ball.p = vector(1, 1, 0)

    i = 0
    while i < n:
        i += 1
        # Randomly choose to move in the x or y direction
        rate(2000)
        while True:
            try:
                dir = random.choice(['x', 'y'])
                pos = random.choice([1, -1])
                if (str(dir) == 'x') and (-side < ball.pos.x + pos < side):
                    ball.pos += vector(pos, 0, 0)
                    break
                elif (str(dir) == 'y') and (-side < ball.pos.y + pos < side):
                    ball.pos += vector(0, pos, 0)
                    break
                else:
                    raise ValueError()
            except ValueError:
                continue

        # ball.pos += ball.p * dt
        # # ball.pos = ball.pos + (ball.p / ball.mass) * dt
        # if not (side > ball.pos.x > -side):
        #     ball.p.x = -ball.p.x
        # if not (side > ball.pos.y > -side):
        #     ball.p.y = -ball.p.y
        # # if not (side > ball.pos.z > -side):
        # #     ball.p.z = -ball.p.z
    # return


if __name__ == '__main__':
    # pb01()
    # pb02()
    pb03()
    pass
