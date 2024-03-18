from rich import print

import numpy as np
### This file contains helper functions for some of the algorithms.
### May later be split into several utils-files


# extended euclidean algorithm for two integers
def x_gcd(a, b):
    if a == 0:
        return b, 0, 1
    gcd, x1, y1  = x_gcd(b % a, a)
    
    x = y1 - (b//a) * x1
    y = x1

    return gcd, x, y

# returns -1*f
def negate_poly(f):
    return [-f[i] for i in range(len(f))]

# returns sum of two polynomials.
# optional modulus n
def poly_add(f, g, n=None):
    if n:
        print(f"Not implemented yet")
        return

    len_f = len(f)
    len_g = len(g)

    if len_f > len_g:
        r = np.copy(f)
        for i in range(len_g):
            r[i] += g[i]
        return r
    if len_f < len_g:
        r = np.copy(g)
        for i in range(len_f):
            r[i] += f[i]
        return r
    if len_f==len_g:
        r = np.array([f[i]+g[i] for i in range(len_f)])
        return r

def poly_sub(f,g):
    return poly_add(f, negate_poly(g))

# straightforward multiplication of polynomials f and g mod ideal x^n + 1
# should use e.g. karatsuba multiplication
# params: 
# polynomials f, g
# n is degree of ideal x^n +1
def poly_mult(f, g, n):
    len_f = len(f)
    len_g = len(g)

    len_fg = len_f+len_g
    f_g = np.full(len_fg, 0)

    for i in range(len_f):
        for j in range(len_g):
            f_g[i+j] += f[i]*g[i]

    # trim zeroes
    f_g = np.trim_zeros(f_g, trim='b')
    #TODO reduce mod x^n +1
    return f_g

# reduce polynomial f mod g over rationals
def poly_reduce_Q(f, g):
    f_rev = np.copy(f)[::-1]
    g_rev = np.copy(g)[::-1]
    quot, rem = np.polydiv(f_rev, g_rev)
    return rem

# check if f is invertible mod p
def is_invertible(f, p):
    # if p is 2, we simply do a parity check
    if p == 2:
        return (np.sum(f)%2)==1

test_f = np.array([1,1,1,1])
test_g = np.array([0,0,1,0])

test_fg = poly_mult(test_f, test_g, 5)