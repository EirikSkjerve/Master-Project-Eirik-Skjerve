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
    log(n), power of degree
    
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
    return ntru_solve(f, g)

def sign(m, priv, logn, seed=None):
    '''
    Calls hawksign
    Inputs: 
    Message m
    Private key priv
    Log2 of degree n
    '''

    np.random.seed(seed)
    rng = RngContext(np.random.randint(0, 256, 40, dtype=np.uint8))

    # turn the message into a numpy array for it to be convertable to bytes
    m_vec = np.array([x for x in m], dtype=str)

    # returns a signature of m
    return hawksign(logn, priv, m_vec, rng)

def verify(m, pub, sig, logn):

    # turn the message into a numpy array for it to be convertable to bytes
    m_vec = np.array([x for x in m], dtype=str)

    return hawkverify(logn, pub, m_vec, sig)

if __name__ == "__main__":
    kgen_s = time()
    keypairs = generate_keypairs(num_pairs=1, logn=8, seed=13)
    kgen_e = time()
    print(f"Keypair generated in {kgen_e - kgen_s} s")
    for i in range(len(keypairs)):

        message = "Hei eg heiter Eirik :)"
        B, priv, pub = keypairs[i]
        
        '''
        f, g, F, G = B
        print(f"f: {f}")
        print(f"g: {g}")
        print(f"F: {F}")
        print(f"G: {G}")
        '''
        sig_s = time()
        signature = sign(message, priv, logn=8, seed=1337)
        sig_e = time()
        print(f"Signature generated in {sig_e-sig_s} s")
        ver_s = time()
        verified = verify(message,pub, signature, logn=8)
        ver_e = time()
        print(f"Signature verified: {verified} in {ver_e - ver_s} s")
        '''
        #print(f"Private key: {priv}")
        #print(f"Public key:  {pub}")
        signature = sign(message, priv, logn=8)
        #print(f"Signature for message: {signature}")
        #print("\n")
        '''