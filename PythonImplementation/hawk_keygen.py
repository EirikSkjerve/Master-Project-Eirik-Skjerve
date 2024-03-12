import numpy as np
from numpy import linalg as LA
import secrets
from utils import *
from rich import print

### High level overview of hawk key-generation ###

# Polynomials will be represented as a list of coefficients, ordered from least significant to most significant term.
# e.g. x^3 + x -1 = [-1, 1, 0, 1]

# encodes an integer to binary
def encode_int(x, k):
    return bin(x)[2:].zfill(k)

# decodes binary to integer
# TODO handle case if input is not a string
def decode_int(bits):
    return int(bits, 2)

# step 1: sample coefficients of f, g through Bin(n). Seed shoud be provided
def sample_coefficients(eta, kgseed):
    
    # set a seed for the sampling
    np.random.seed(kgseed)

    # numpy's binomial distribution. Stored as np-array of type int
    # samples are centered with -eta/2
    centred_samples = np.array(np.random.binomial(eta,p=0.5,size=n) - eta/2, dtype=np.int8)
    
    return centred_samples

# generate f and g
def generate_f_g():
    
    # eta for n=256 is 2
    eta = 2

    kgseed_f = 13
    f = sample_coefficients(eta, kgseed_f)
    kgseed_g = 27
    g = sample_coefficients(eta, kgseed_g)

    print(f"f: {f}\ng:{g}")
    # restart with a new seed if conditions are not fulfilled
    return (f, g)

# step 2: check conditions for f, g
def verify_f_g(f, g) -> bool:
    
    norm_cond = 20  # some threshold for the norms of f, g
    f_norm = LA.norm(f)
    g_norm = LA.norm(g)
    
    if f_norm > norm_cond or g_norm > norm_cond:
        print(f"Norm too big!")
        return False

    #TODO check some other properties aswell

    return True

# step 3: get r = NTRUSolve(f, g) = F, G s.t. fG -gF = 1 mod ideal. This is the NTRU-equation
# see https://github.com/hawk-sign/hawk-py/blob/main/ntrugen/ntrugen_hawk.py
def NTRU_solve(f, g):

    n = len(f)
    if n == 1:
        f0 = f[0]
        g0 = g[0]
        d, u, v = x_gcd(f0, g0)
        if d != 1:
            raise ValueError
        else:
            return [-q * v], [q * u]

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

# step 9: return (pk, sk) = (Q, (B, hpub)). Guessing B is the private key.

def hawk_keygen():
    n = 256
    ideal = np.array([1] + ([0]* (n-2)) + [1], dtype=np.uint8)

    kgseed = np.random()
    f, g = generate_f_g(kgseed)

if __name__ == "__main__":
    # test environment

    # degree n
    n = 256
    q = 1
    ideal = np.array([1] + ([0]* (n-2)) + [1], dtype=np.uint8)
    print(f"Ideal: {ideal}")

    # for hawk256, Bin(eta) param is 2
    eta = 2
    # use 'secrets' module for CSPRNG
    #kgseed = secrets.randbits(5)