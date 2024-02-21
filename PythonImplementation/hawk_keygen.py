import numpy as np
from Crypto.Hash import SHAKE256
from typing import *

### High level overview of hawk key-generation ###

# step 1: sample coefficients of f, g through Bin(n)
def sample_coefficients(n) -> tuple:
    pass

# step 2: check conditions for f, g
def f_g_conditions(f, g) -> bool:
    pass

# step 3: get r = NTRUSolve(f, g)
def NTRU_solve(f, g):
    pass

# step 4: check orthogonality. If r is orthogonal, restart
def check_orth(r) -> bool:
    pass

# step 5: set (F, G) = r

# step 6: set B = Matrix 
#(f, F)
#(g, G),
# set Q = B*B

# step 7: check KGen-encoding of Q and B. Restart if condition fails

# step 8: set hpub to H(Q) (hash function)

# step 9: return (pk, sk) = (Q, (B, hpub)). Guessing Q is the private key and (B, hpub) is the public key.