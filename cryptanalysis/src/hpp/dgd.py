import math

s = 2.02


def rho(x):
    return math.exp((-x**2)/(2*s**2))


def pr(x):
    return rho(x)/sum([rho(y) for y in range(-100, 100)])


def expectance():
    return sum([pr(x)*(x**2) for x in range(-100, 100)])


def fourth():
    return sum([pr(x)*(x**4) for x in range(-100, 100)])


if __name__ == "__main__":
    print(f"s = {s}")
    print(f"s^2 = {s**2}")
    print(f"E[x^2] = {expectance()}")
    print(f"E[x^4] = {fourth()}")
    print(f"3*s^4 = {3*s**4}")
