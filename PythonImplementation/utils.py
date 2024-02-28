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
