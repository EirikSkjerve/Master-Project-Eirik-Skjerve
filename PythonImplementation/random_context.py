import numpy as np

import random

class RandomContext:
    '''
    Generates random numbers in a deterministic way
    '''
    # initialize context with a base seed
    def __init__(self, base_seed):
        self.base_seed = base_seed
        self.i = 0

    def random(self, bits):
        random.seed((self.base_seed + self.i))
        self.i += 1
        randbits = random.getrandbits(bits)
        return randbits
