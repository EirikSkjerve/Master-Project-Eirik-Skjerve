# we import HAWK's python implementation of the hawk key generation algorithm
from hawk_py.keygen import hawkkeygen, RngContext
from hawk_py.ntrugen.ntrugen_hawk import ntru_solve
import numpy as np

def generate_keypairs(num_pairs, logn=8, seed=None):
    '''
    Input: 
    Wanted number of pairs
    log(n), power of degree
    
    Outut:
    Randomly generated public/private key pairs for HAWK
    '''
    
    # TODO set an RNGcontext object from seed for testing

    #rng = RngContext(seed)

    keypairs = [None]*num_pairs
    for p in range(num_pairs):
        np.random.seed(seed+p)
        rng = RngContext(np.random.randint(0, 256, 40, dtype=np.uint8))
        keypairs[p] = hawkkeygen(logn, rng)
        print(f"{p+1}/{num_pairs} generated")
    return keypairs

def generate_F_G(f, g):
    return ntru_solve(f, g)

if __name__ == "__main__":
    keypairs = generate_keypairs(2, seed=13)
    for i in range(len(keypairs)):

        priv, pub = keypairs[i]
        print(f"Private key: {priv}")
        print(f"Public key:  {pub}")