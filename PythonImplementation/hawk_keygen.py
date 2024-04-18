import numpy as np
from numpy import linalg as LA
from utils import *
from rich import print
from random_context import RandomContext
import random
# currently using this as a black-box
from NTRUSolve import NTRU_solve_bruteforce
### High level overview of hawk key-generation ###

# Polynomials will be represented as a list of coefficients, ordered from least significant to most significant term.
# e.g. x^3 + x -1 = [-1, 1, 0, 1]

# encodes an integer to binary representation in a numpy array
def encode_int(x, k):
    encoded = np.full(k, 0, dtype=np.uint8)
    bin_string = bin(x)[2:].zfill(k)
    for i, b_s in enumerate(bin_string):
        encoded[i] = int(b_s)
    return encoded

# decodes binary to integer
def decode_int(bits):
    return int(''.join([str(b) for b in bits]), 2)

# step 1: sample coefficients of f, g through Bin(n). Seed shoud be provided
def sample_coefficients(eta, coef_seed, n):
    
    # set a seed for the sampling
    np.random.seed(coef_seed)

    # numpy's binomial distribution. Stored as np-array of type int
    # samples are centered with -eta/2
    centred_samples = np.array(np.random.binomial(2*eta,p=0.5,size=n) - eta, dtype=np.int8)
    
    return centred_samples

# generate f and g
def generate_f_g(random_context, n):
    
    # eta for n=256 is 2
    eta = 2

    f = sample_coefficients(eta, random_context.random(10), n)

    g = sample_coefficients(eta, random_context.random(10), n)

    print(f"f: {f}\ng: {g}")

    return (f, g)

# step 2: check conditions for f, g
def verify_f_g(f, g) -> bool:
    
    threshold = 256*(1.042**2)
    f_norm = LA.norm(f)
    g_norm = LA.norm(g)
    
    if f_norm > threshold or g_norm > threshold: 
        print(f"Norm too big!")
        return False

    return True

# step 5: set (F, G) = r

# step 6: set B = Matrix 
#(f, F)
#(g, G),
# set Q = B*B

# step 7: check KGen-encoding of Q and B. Restart if condition fails

# step 8: set hpub to H(Q) (hash function)

# step 9: return (pk, sk) = (Q, (B, hpub)). Guessing B is the private key.

def hawk_keygen(retry=False):
    n = 64

    global ideal
    ideal = np.array([1] + ([0]* (n-2)) + [1], dtype=np.uint8)

    base_seed = 14

    if retry:
        base_seed = random.randint(100)

    random_context = RandomContext(base_seed)

    f, g = generate_f_g(random_context, n)

    print(f"f: {list(f)}")
    print(f"g: {list(g)}")
    # retry if f and g does not meet requirements
    if not verify_f_g(f, g):
        return hawk_keygen(retry=True)

    # trying to solve ntru-equation
    F, G = NTRU_solve_bruteforce(f, g)
    #F, G = NTRU_solve(list(f), list(g))

    quit()
    assert (int(poly_reduce_Q(poly_mult(f, G, ideal) - poly_mult(g,F, ideal))) == 1)

    # B is the secret basis for the lattice
    B = np.array([[f, F],
                 [g, G]])

    # Inverse of B
    B_inv = np.array([[G, negate_poly(F)],
                      [negate_poly(g)], f])

    # Q is the public key. Excplicitly it is computed the following way: 
    Q = [[None, None],
         [None, None]]

    Q[0][0] = poly_add(poly_mult(f, f), poly_mult(g, g))
    Q[0][1] = poly_add(poly_mult(f, F), poly_mult(g, G))
    Q[1][0] = poly_add(poly_mult(F, f), poly_mult(G, g))
    Q[1][1] = poly_add(poly_mult(F, F), poly_mult(G, G))

    return

if __name__ == "__main__":
    # test environment
    pub, priv = hawk_keygen()