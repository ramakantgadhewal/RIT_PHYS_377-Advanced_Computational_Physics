"""
Author: Ramsey (Rayla) Phuc
"""
import random

import matplotlib.pyplot as plt
import numpy as np


def pb01():
    # Generate a (n x n) matrix to represent a system of (n x n) spin on a square lattice
    n = 20
    grid = np.zeros((n, n), dtype=float)

    # Randomly set each grid value to be either 1 or -1
    choices = [-1, 1]
    for r in range(0, n):
        for c in range(0, n):
            grid[r][c] = random.choice(choices)

    """
    Part (a): First write a function to calculate the total energy of the system, as 
    given by the equation above. That is, for a given array of values of the spins, go 
    through every pair of adjacent spins and add up the contributions s_i s_j from all of 
    them, then multiply by −J.
    """

    def getEnergy(grid):
        """
        Calculates the total energy of the system.
        """
        return np.sum(grid[:, 0:n - 1] * grid[:, 1:n]) + \
               np.sum(grid[0:n - 1, :] * grid[1:n, :])

    J = 1

    """
    Part (b): Now use your function as the basis for a Metropolis-style simulation of 
    the Ising model with J = 1 and temperature T = 1 in units where the Boltzmann 
    constant kB is also 1. Initially set the spin variables randomly to ±1, so that 
    on average about a half of them are up and a half down, giving a total 
    magnetization of roughly zero. Then choose a spin at random, flip it, and calculate 
    the new energy after it is flipped, and hence also the change in energy as a result 
    of the flip. Then decide whether to accept the flip using the Metropolis acceptance 
    formula. If the move is rejected you will have to flip the spin back to where it 
    was. Otherwise you keep the flipped spin. Now repeat this process for many moves.
    """
    # Constants and Parameters
    T = 1  # Temperature
    kb = 1  # Boltzmann Constant
    steps = 1_000_000

    # Calculate the current energy of the system
    Ei = -J * getEnergy(grid)

    # Store the total energy
    E_total = [Ei]

    # Record the Magnetization values
    M = [np.sum(grid)]

    # Repeat steps times
    for q in range(0, steps):
        # Pick a random point on the grid
        r, c = random.randint(0, n - 1), random.randint(0, n - 1)
        # Flip the spin
        grid[r][c] *= -1
        # Calculate the Energy of this grid
        Ej = -J * getEnergy(grid)
        # Calculate the change in energy
        dE = Ej - Ei
        # Decide to accept or reject the flip
        # If accept: Then set the old energy equal to the current energy
        # Otherwise, re-flip the spin
        if dE <= 0:
            Ei = Ej
        else:
            if random.random() > np.exp(-dE / (kb * T)):
                grid[r][c] *= -1
            else:
                Ei = Ej
        # Record current energy and magnetization values
        E_total.append(Ei)
        M.append(np.sum(grid))
    """
    Part (c): Make a plot of the total magnetization of the system as a function of 
    time for a million Monte Carlo steps. You should see that the system develops a 
    “spontaneous magnetization,” a nonzero value of the overall magnetization.
    """
    # Plotting the Magnetization and the Energy
    fig = plt.figure()
    plt.plot(E_total, color="black", label="Energy")
    plt.plot(M, color="red", label="Magnetization")
    plt.minorticks_on()
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    plt.legend()
    plt.xlabel("Number of Steps")
    plt.ylabel("Energy/Magnetization")
    fig.savefig('hw09.png')
    plt.show()

    """
    Part (d): Run your program several times and observe the sign of the magnetization 
    that develops, positive or negative. Describe what you find and give a brief 
    explanation of what is happening.
    
    Answer: Running the program several times would yield drastically different results. 
    Sometimes the magnetization would rapidly converge to -400 and other times it would 
    rapidly converge to +400. 
    """
    # Plotting the model to verify that the above code works. The goal of this chunk of
    # code is to make sure if we have a plot of all (or mostly one color). If this is
    # true, then this implies that all (if not most) of the spins are in one direction.

    # x, y = np.meshgrid(np.arange(0, n, 1), np.arange(0, n, 1))
    # plt.contourf(x, y, grid)
    # plt.show()
    return


if __name__ == '__main__':
    pb01()
    pass
