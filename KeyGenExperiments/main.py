# we import HAWK's python implementation of the hawk key generation algorithm
from hawk_py.keygen import hawkkeygen, RngContext
from hawk_py.sign import hawksign
from hawk_py.ntrugen.ntrugen_hawk import ntru_solve
from hawk_py.verify import hawkverify

from time import time
import numpy as np

def generate_keypairs(num_pairs, logn, seed=None):
    '''
    Input: 
    Wanted number of pairs
    log2 of degree n
    seed (optional)
    
    Outut:
    Randomly generated public/private key pairs for HAWK
    '''
    
    # initialize empty list for keypairs
    keypairs = [None]*num_pairs
    for p in range(num_pairs):
        # set seed
        np.random.seed(seed+p)

        # create an RngContext used for hawkkeygen
        rng = RngContext(np.random.randint(0, 256, 40, dtype=np.uint8))

        # create a keypair from HAWK's implementation
        keypairs[p] = hawkkeygen(logn, rng)
        print(f"{p+1}/{num_pairs} generated")

    return keypairs

def generate_F_G(f, g):
    '''
    Calls ntru_solve on f and g
    Inputs:
    polynomials f, g

    Outputs:
    polynomials F, G such that fG - gF = 1 mod (x^n +1)
    '''
    return ntru_solve(f, g)

def sign(m, priv, logn, seed=None):
    '''
    Calls hawksign
    Inputs: 
    Message m
    Private key priv
    Log2 of degree n
    
    Outputs:
    Signature sig
    '''

    np.random.seed(seed)
    rng = RngContext(np.random.randint(0, 256, 40, dtype=np.uint8))

    # turn the message into a numpy array for it to be convertable to bytes
    m_vec = np.array([x for x in m], dtype=str)

    # returns a signature of m
    return hawksign(logn, priv, m_vec, rng)

def verify(m, pub, sig, logn):
    '''
    Calls hawkverify
    Inputs:
    Message m
    Public key pub
    Signature sig
    Log2 of degree n

    Outputs: 
    True or False
    '''

    # turn the message into a numpy array for it to be convertable to bytes
    m_vec = np.array([x for x in m], dtype=str)

    return hawkverify(logn, pub, m_vec, sig)

if __name__ == "__main__":
    # generate keypair(s)

    num_pairs = 3
    logn = 10

    kgen_s = time()
    keypairs = generate_keypairs(num_pairs, logn, seed=13)
    kgen_e = time()
    print(f"{num_pairs} keypair(s) generated in {kgen_e - kgen_s} s")

    # run a basic pipeline
    for i in range(len(keypairs)):

        message = "A test message for HAWK digital signature scheme"
        B, priv, pub = keypairs[i]
        
        f, g, F, G = B
        # print f, g, F, G
        '''
        print(f"f: {f}")
        print(f"g: {g}")
        print(f"F: {F}")
        print(f"G: {G}")
        '''
        
        # signature generation
        sig_s = time()
        signature = sign(message, priv, logn, seed=1337)
        sig_e = time()
        print(f"Signature generated in {sig_e-sig_s} s")

        # signature verification
        ver_s = time()
        verified = verify(message,pub, signature, logn)
        ver_e = time()
        print(f"Signature verified: {verified} in {ver_e - ver_s} s")
        print("\n")