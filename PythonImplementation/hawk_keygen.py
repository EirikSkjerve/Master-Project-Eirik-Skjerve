import numpy as np
#from Crypto.Hash import SHAKE256
import secrets

### High level overview of hawk key-generation ###
def encode_int(x, k):
    bin_rep = bin(x)[2:].zfill(k)
    print(f"{int(bin_rep, 2)}")
    return bin_rep

def decode_int(bits, k):
    return int(bits, 2)

# step 1: sample coefficients of f, g through Bin(n)
def sample_coefficients(eta, kseed):
    np.random.seed(kseed)
    centred_samples = np.random.binomial(eta,p=0.5,size=n) - eta/2
    
    return centred_samples

def generate_f_g(eta):
    kseed_f = 13
    f = sample_coefficients(eta, kseed_f)
    kseed_g = 27
    g = sample_coefficients(eta, kseed_g)

    # restart with a new seed if conditions are not fulfilled
    return (f, g) if verify_f_g else generate_f_g(eta)

# step 2: check conditions for f, g
def verify_f_g(f, g) -> bool:
    return True

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

# step 9: return (pk, sk) = (Q, (B, hpub)). Guessing B is the private key.


if __name__ == "__main__":
    # test environment

    # degree n
    n = 256

    # for hawk256, Bin(eta) param is 2
    eta = 2
    # use 'secrets' module for CSPRNG
    #kseed = secrets.randbits(5)
    kseed = 13
    sample_coefficients(eta, kseed)