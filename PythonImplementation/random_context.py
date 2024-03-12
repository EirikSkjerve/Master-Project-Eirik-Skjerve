import numpy as np
import secrets

class RandomContext:
    # initialize context with a base seed
    def __init__(self, base_seed):
        self.base_seed = base_seed