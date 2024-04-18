# we import HAWK's python implementation of the hawk key generation algorithm
from hawk_py.keygen import hawkkeygen
from hawk_py.ntrugen.ntrugen_hawk import ntru_solve
def generate_keypairs(num_pairs, logn=8):
    '''
    Input: 
    Wanted number of pairs
    log(n), power of degree
    
    Outut:
    Randomly generated public/private key pairs for HAWK
    '''
    keypairs = [None]*num_pairs
    for p in range(num_pairs):
        keypairs[p] = hawkkeygen(logn)
        print(f"{p+1}/{num_pairs} generated")
    return keypairs

def generate_F_G(f, g):
    return ntru_solve(f, g)

if __name__ == "__main__":
    keypairs = generate_keypairs(10)
    print(keypairs)