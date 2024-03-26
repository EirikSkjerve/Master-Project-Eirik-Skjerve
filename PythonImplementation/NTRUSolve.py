from utils import x_gcd, poly_mult, poly_reduce_Q

# these three functions are borrowed from ntrugen_hawk.py
# currently not using the fft or karatsuba for multiplication
def lift(a):
    '''
    Lift an element a of Q[x] / (x ** (n//2)+1) up to Q[x] / (x**n + 1)
    '''
    n = len(a)
    lifted = [0]*(2*n)
    for i in range(n):
        lifted[2*i] = a[i]
    return lifted

def galois_conjugate(a):
    """
    Galois conjugate of an element a in Q[x] / (x ** n + 1).
    Here, the Galois conjugate of a(x) is simply a(-x)
    """
    n = len(a)
    # inserted %2 to reduce calculations
    return [((-1) ** (i%2)) * a[i] for i in range(n)]


def field_norm(a):
    """
    Project an element a of Q[x] / (x ** n + 1) onto Q[x] / (x ** (n // 2) + 1).
    Only works if n is a power-of-two.
    """
    n2 = len(a) // 2
    ideal = [1]+[0]*(len(a)-2)+[1]
    ae = [a[2 * i] for i in range(n2)]
    ao = [a[2 * i + 1] for i in range(n2)]

    ae_squared = poly_mult(ae, ae, ideal)
    ao_squared = poly_mult(ao, ao, ideal)
    res = ae_squared[:]

    for i in range(n2 - 1):
        res[i + 1] -= ao_squared[i]
    res[0] += ao_squared[n2 - 1]
    return res

# step 3: get r = NTRUSolve(f, g) = F, G s.t. fG -gF = 1 mod ideal. This is the NTRU-equation
# see https://github.com/hawk-sign/hawk-py/blob/main/ntrugen/ntrugen_hawk.py
def NTRU_solve(f, g):

    q = 1
    
    n = len(f)
    ideal = [1]+[0]*(n-2)+[1]
    if n == 1:
        f0 = f[0]
        g0 = g[0]
        d, u, v = x_gcd(f0, g0)
        if d != 1:
            raise ValueError
        else:
            return [-q * v], [q * u]
    else:
        fp = field_norm(f)
        gp = field_norm(g)
        Fp, Gp = NTRU_solve(fp, gp)
        F = poly_mult(lift(Fp), galois_conjugate(g), ideal)
        G = poly_mult(lift(Gp), galois_conjugate(f), ideal)
        return F, G


