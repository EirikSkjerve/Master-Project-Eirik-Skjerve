import math
import matplotlib.pyplot as plt


def rho(x, s):
    return math.exp((-x**2)/(2*s**2))


def pr(x, s, eta, c):
    return rho(x, s)/sum([rho(y, s) for y in range(-eta+c, eta+c, 2)])


def pdf(s, eta, c):
    return sum([pr(x, s, eta, c) for x in range(-eta+c, eta+c, 2)])


# def pdf2(s, eta):
#     return sum([rho(x, s) / sum([rho(y, s) for y in range(-eta, eta)]) for x in range(-eta, eta)])
#

def expectance(s, eta, c):
    return sum([pr(x, s, eta, c)*(x**2) for x in range(-eta+c, eta+c, 2)])


def fourth(s, eta, c):
    return sum([pr(x, s, eta, c)*(x**4) for x in range(-eta+c, eta+c, 2)])


if __name__ == "__main__":

    # sigma
    s = 1.01
    eta_range = 20
    # c = 0

    # diffs_1 = []
    # diffs_2 = []
    # for i in range(0, eta_range):
    #     diffs_1.append(abs(expectance(s, i)-s**2))
    #     diffs_2.append(abs(fourth(s, i)-3*s**4))
    #
    # plt.grid(True)
    # plt.plot(range(0, eta_range), diffs_1, color="|green", label="Exp(x^2)")
    # plt.plot(range(0, eta_range), diffs_2, label="|Exp(x^4)")
    #
    # plt.xlabel("eta")
    # plt.ylabel("2nd and 4th moment distance from s^2 and 3s^4 respectively")
    #
    # plt.legend()
    #
    # plt.show()

    print(f"\ns = {s}")
    t1 = 0
    t2 = 0
    for c in range(2):

        print(f"\ns^2 = {s**2}")
        print(f"E[x^2] = {expectance(s, eta_range, c)}")
        # diff_1 = abs(expectance(s, eta_range, c) - s**2)
        t1 += expectance(s, eta_range, c)
        # print(f"diff = {diff_1} \n")
        print(f"3*s^4 = {3*s**4}")
        print(f"E[x^4] = {fourth(s, eta_range, c)} \n")
        # diff_2 = abs(fourth(s, eta_range, c) - 3*s**4)
        t2 += fourth(s, eta_range, c)
        # print(f"diff = {diff_2}")

    print(f"Avg. E[x^2]: {t1/2}")
    print(f"Avg. E[x^4]: {t2/2}")
    print(f"3E[x^2]^2: {3*(t1/2)**2}")
