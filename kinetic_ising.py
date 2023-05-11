#!/usr/bin/env python3
"""
GPLv3 2020 Miguel Aguilera

This code allows to run simulations of the kinetic Ising model,
with asymmetric weights and parallel updates.
"""

import numpy as np

def bool2int(x):          # Transform bool array into positive integer
    y = 0
    for i, j in enumerate(np.array(x)[::-1]):
        y += j * 2**i
    return y


def bitfield(n, size):  # Transform positive integer into bit array
    x = [int(x) for x in bin(n)[2:]]
    x = [0] * (size - len(x)) + x
    return np.array(x)


class ising:                        # Asymmetric Ising model simulation class

    def __init__(self, netsize):  # Create ising model

        self.size = netsize                        # Network size
        self.H = np.zeros(netsize)                # Fields
        self.J = np.zeros((netsize, netsize))  # Couplings
        self.Beta = 1                            # Inverse temperature

        self.randomize_state()                    # Set random state

    def randomize_state(self):        # Randomize network state
        self.s = np.random.randint(0, 2, self.size) * 2 - 1

    def random_fields(self):        # Set random values for H
        self.H = np.random.rand(self.size) * 2 - 1

    def random_wiring(self):        # Set random values for J
        self.J = np.random.randn(self.size, self.size) / self.size

    # Update the state of the network using Little parallel update rule
    def ParallelUpdate(self):
        self.h = self.H + np.dot(self.J, self.s)
        r = np.random.rand(self.size)
        self.s = -1 + 2 * (2 * self.Beta * self.h > -
                           np.log(1 / r - 1)).astype(int)
                           
    def GlauberStep(self):                  # Execute step of the Glauber algorithm
        i = np.random.randint(self.size)
        h = self.H[i] + np.dot(self.J[i,:],self.s)
        self.s[i] = int(np.random.rand()*2-1 < np.tanh(h))*2-1   # Glauber

    def SequentialGlauberStep(self):        # Execute N steps of the Glauber algorithm
        for i in range(self.size):
            self.GlauberStep()

