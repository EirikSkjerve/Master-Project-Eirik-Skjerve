Hawk256_T0 = [
    0x26B871FBD58485D45050,
    0x07C054114F1DC2FA7AC9,
    0x00A242F74ADDA0B5AE61,
    0x0005252E2152AB5D758B,
    0x00000FDE62196C1718FC,
    0x000000127325DDF8CEBA,
    0x0000000008100822C548,
    0x00000000000152A6E9AE,
    0x0000000000000014DA4A,
    0x0000000000000000007B,
]

Hawk512_T0 = [
    0x2C058C27920A04F8F267,
    0x0E9A1C4FF17C204AA058,
    0x02DBDE63263BE0098FFD,
    0x005156AEDFB0876A3BD8,
    0x0005061E21D588CC61CC,
    0x00002BA568D92EEC18E7,
    0x000000CF0F8687D3B009,
    0x0000000216A0C344EB45,
    0x0000000002EDF0B98A84,
    0x0000000000023AF3B2E7,
    0x00000000000000EBCC6A,
    0x000000000000000034CF,
    0x00000000000000000006,
]

Hawk1024_T0 = [
    0x2C583AAA2EB76504E560,
    0x0F1D70E1C03E49BB683E,
    0x031955CDA662EF2D1C48,
    0x005E31E874B355421BB7,
    0x000657C0676C029895A7,
    0x00003D4D67696E51F820,
    0x0000014A1A8A93F20738,
    0x00000003DAF47E8DFB21,
    0x0000000006634617B3FF,
    0x000000000005DBEFB646,
    0x00000000000002F93038,
    0x0000000000000000D5A7,
    0x00000000000000000021,
]

Hawk256_T1 = [
    0x13459408A4B181C718B1,
    0x027D614569CC54722DC9,
    0x0020951C5CDCBAFF49A3,
    0x0000A3460C30AC398322,
    0x000001355A8330C44097,
    0x00000000DC8DE401FD12,
    0x00000000003B0FFB28F0,
    0x00000000000005EFCD99,
    0x00000000000000003953,
    0x00000000000000000000,
]

Hawk512_T1 = [
    0x1AFCBC689D9213449DC9,
    0x06EBFB908C81FCE3524F,
    0x01064EBEFD8FF4F07378,
    0x0015C628BC6B23887196,
    0x0000FF769211F07B326F,
    0x00000668F461693DFF8F,
    0x0000001670DB65964485,
    0x000000002AB6E11C2552,
    0x00000000002C253C7E81,
    0x00000000000018C14ABF,
    0x0000000000000007876E,
    0x0000000000000000013D,
    0x00000000000000000000,
]

Hawk1024_T1 = [
    0x1B7F01AE2B17728DF2DE,
    0x07506A00B82C69624C93,
    0x01252685DB30348656A4,
    0x001A430192770E205503,
    0x00015353BD4091AA96DB,
    0x000009915A53D8667BEE,
    0x00000026670030160D5F,
    0x00000000557CD1C5F797,
    0x00000000006965E15B13,
    0x00000000000047E9AB38,
    0x000000000000001B2445,
    0x000000000000000005AA,
    0x00000000000000000000,
]

Hawk256_T0_unscaled = [x / (2**78) for x in Hawk256_T0]
Hawk256_T1_unscaled = [x / (2**78) for x in Hawk256_T1]

Hawk512_T0_unscaled = [x / (2**78) for x in Hawk512_T0]
Hawk512_T1_unscaled = [x / (2**78) for x in Hawk512_T1]

Hawk1024_T0_unscaled = [x / (2**78) for x in Hawk1024_T0]
Hawk1024_T1_unscaled = [x / (2**78) for x in Hawk1024_T1]


import random
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
from decimal import Decimal, getcontext
import math
from mpl_toolkits.mplot3d import Axes3D

def rho(x, sigma):
    return math.e ** (-(x**2) / (2 * (sigma**2)))


def pr(x, sigma):
    return (rho(x, sigma)) / sum([rho(z, sigma) for z in range(-100, 100)])


def sample(c):
    q = random.uniform(-1, 1)

    T = None
    if c == 0:
        T = Hawk256_T0_unscaled
    if c == 1:
        T = Hawk256_T1_unscaled

    z = 0
    v = 0
    while T[z] != 0:
        if abs(q) <= T[z]:
            v += 1
        z += 1
        if z == len(T):
            break

    v = 2 * v + c
    if q < 0:
        v = -v
    return v


def sample_vector(t):
    x = []
    for c in t:
        x.append(sample(c))

    return x


def barchart(x):
    # Count occurrences
    freqs = Counter(x)
    num_samples = len(x)

    # Extract values for plotting
    x_values = list(freqs.keys())
    y_values = [x / num_samples for x in list(freqs.values())]
    print(sum(y_values))

    # Plot the bar chart
    plt.bar(x_values, y_values)

    # Labels and title
    plt.xlabel("Values")
    plt.ylabel("Relative frequencies")
    plt.title("Bar chart for the discrete gaussian distribution with sigma = 2.02")

    # Ensure all x-tick labels are visible
    plt.xticks(x_values)  # Explicitly set the x-ticks
    # Show the plot
    plt.show()


def scatterplot():
    num_samples = 100000
    V = np.array([[3, 1], [2, 2]])

    l = np.linalg.cholesky(np.linalg.inv(V.T @ V))
    C = l.T @ V.T
    # C = np.array([[1, 0],
    #              [0, 1]])
    print(V)
    print(C)
    ws = []
    cs = []
    for _ in range(num_samples):
        t = [random.randint(0, 1), random.randint(0, 1)]
        x = sample_vector(t)
        ws.append(V.T @ x)
        cs.append(C.T @ x)

    points = np.array(cs)
    x = points[:, 0]
    y = points[:, 1]

    v1 = C[0, :]  # First column
    v2 = C[1, :]  # Second column

    # u = np.array([0.8944, 0.4472])
    # u /= np.linalg.norm(u)
    # print(np.mean(((points / 2.02) @ u)**4))

    plt.scatter(x, y, color="black", s=0.1, label='transformed signature sample')

    plt.arrow(
        0,
        0,
        v1[0],
        v1[1],
        head_width=0.1,
        width=0.05,
        length_includes_head=True,
        color="red",
    )

    plt.arrow(
        0,
        0,
        v2[0],
        v2[1],
        head_width=0.1,
        width=0.05,
        length_includes_head=True,
        color="blue",
    )

    # Optionally, add labels and title
    plt.xlabel("X")
    plt.ylabel("Y")
    # Legg til et kjapt lite “hack” for å vise ‘label’ i legend
    # ved å “tegne” en usynlig linje og kalle den 'label'
    plt.plot([], [], "r-", label="C column 1")
    plt.plot([], [], "b-", label="C column 2")
    plt.legend()

    plt.axis('equal')

    plt.grid(True)
    # Show the plot
    plt.show()


def scatterplot_3d():
    num_samples = 10000

    # Example 3x3 matrix (adapt as needed)
    V = np.array([
        [3, 1, 1],
        [2, 2, 1],
        [1, 1, 3]
    ])

    # Cholesky-based transformation, analogous to your 2D code
    L = np.linalg.cholesky(np.linalg.inv(V.T @ V))
    C = L.T @ V.T

    # Store transformed points
    cs = []

    for _ in range(num_samples):
        # Example: pick a random [0 or 1] for each dimension
        # Replace this with your own 3D sample_vector() if preferred
        t = [random.randint(0, 1),
             random.randint(0, 1),
             random.randint(0, 1)]
        x = sample_vector(t)

        cs.append(C.T @ x)  # 3D transformed point

    points = np.array(cs)
    x_vals = points[:, 0]
    y_vals = points[:, 1]
    z_vals = points[:, 2]

    # Extract vectors from C
    v1 = C[0, :]
    v2 = C[1, :]
    v3 = C[2, :]

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the scatter of points
    ax.scatter(x_vals, y_vals, z_vals, color="black", s=0.1)

    # Draw the three vectors as arrows
    ax.quiver(
        0, 0, 0, 
        v1[0], v1[1], v1[2],
        length=1, color='red'
    )
    ax.quiver(
        0, 0, 0, 
        v2[0], v2[1], v2[2],
        length=1, color='blue'
    )
    ax.quiver(
        0, 0, 0, 
        v3[0], v3[1], v3[2],
        length=1, color='green'
    )

    # Label axes
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    plt.title("3D Scatter Plot Example")
    plt.show()


getcontext().prec = 1000

T0_256 = [Decimal(x) / Decimal(2**78) for x in Hawk256_T0]
T1_256 = [Decimal(x) / Decimal(2**78) for x in Hawk256_T1]

T0_512 = [Decimal(x) / Decimal(2**78) for x in Hawk512_T0]
T1_512 = [Decimal(x) / Decimal(2**78) for x in Hawk512_T1]

T0_1024 = [Decimal(x) / Decimal(2**78) for x in Hawk1024_T0]
T1_1024 = [Decimal(x) / Decimal(2**78) for x in Hawk1024_T1]


# Pr(X = x) in DGD
def pdf(x, n):
    c = x % 2

    match n:
        case 256:
            T0, T1 = T0_256, T1_256
        case 512:
            T0, T1 = T0_512, T1_512
        case 1024:
            T0, T1 = T0_1024, T1_1024

    match c:
        case 0:
            T = T0
        case 1:
            T = T1

    if abs(x) == 0:
        return Decimal(1 / 2) * Decimal(1 - T[0])

    if abs(x) == 1:
        return Decimal(1 / 4) * Decimal(1 - T[0])

    if abs(x) - c == len(T) * 2:
        return Decimal(1 / 4) * Decimal(T[round(((abs(x) - c) / 2)) - 1])

    return Decimal(1 / 4) * Decimal(
        T[round(((abs(x) - c) / 2)) - 1] - T[round((abs(x) - c) / 2)]
    )


def expectation_x(n):
    match n:
        case 256:
            interval = 2 * 10 + 1
        case 512:
            interval = 2 * 13 + 1
        case 1024:
            interval = 2 * 13 + 1
    res = Decimal(0)
    for x in range(-interval, interval + 1):
        temp = Decimal(x) * pdf(x, n)
        res += temp
    return res


def variance_x(n):
    match n:
        case 256:
            interval = 2 * 10 + 1
        case 512:
            interval = 2 * 13 + 1
        case 1024:
            interval = 2 * 13 + 1
    res = Decimal(0)
    for x in range(-interval, interval + 1):
        temp = Decimal(x**2) * pdf(x, n)
        res += temp
    return res


def sigma_x(n):
    return variance_x(n) ** Decimal(0.5)


def kurtosis_x(n):
    match n:
        case 256:
            interval = 2 * 10 + 1
        case 512:
            interval = 2 * 13 + 1
        case 1024:
            interval = 2 * 13 + 1
    res = Decimal(0)
    for x in range(-interval, interval + 1):
        temp = Decimal(x**4) * pdf(x, n)
        res += temp
    return res


def expectation_z(n):
    return Decimal(0)


def variance_z(n):
    return Decimal(1)


def sigma_z(n):
    return variance_z(n) ** Decimal(0.5)


def kurtosis_z(n, sig4):
    return kurtosis_x(n) / sig4


def mom_k_x(n, k):
    match n:
        case 256:
            interval = 2 * 10 + 1
        case 512:
            interval = 2 * 13 + 1
        case 1024:
            interval = 2 * 13 + 1
    res = Decimal(0)
    for x in range(-interval, interval + 1):
        temp = Decimal(x**k) * pdf(x, n)
        res += temp
    return res


def mom_k_z(n, k, sig_x):
    return mom_k_x(n, k) / sig_x ** Decimal(k)


def norm_mom_k(k):
    import math
    if k % 2 == 1:
        return 0

    return (2**(-k/2)) * (math.factorial(k) / (math.factorial(round(k/2))))


def pr_y_only(k, n):
    return pdf(2*k, n) + pdf((2*k) + 1, n)


def exp_y(n):
    exp = 0
    for i in range(-10, 11):
        exp += i*pr_y_only(i, 256)

    return exp


def var_y(n):
    var = 0
    for i in range(-10, 11):
        var += (i**2)*pr_y_only(i, 256)

    return var


def sig_y(n):
    return math.sqrt(var_y(n))


def kur_y(n):
    kur = 0
    for i in range(-10, 11):
        kur += (i**4)*pr_y_only(i, 256)

    return kur

if __name__ == "__main__":

    n = 256

    exp_x = expectation_x(n)
    sigma_x = sigma_x(n)
    var_x = variance_x(n)
    kur_x = kurtosis_x(n)
    exp_z = expectation_z(n)
    sigma_z = sigma_z(n)
    var_z = variance_z(n)
    kur_z = kurtosis_z(n, sigma_x ** Decimal(4))

    random_mom = mom_k_z(n, 20, sigma_x)
    # print(random_mom)
    # print(norm_mom_k(20))

    sum = 0
    for i in range(-10, 11):
        sum += pr_y_only(i, n)

    print(sum)

    expy = exp_y(n)
    vary = var_y(n)
    sigy = sig_y(n)
    kury = kur_y(n)

    print(f"E[y] = {expy}")
    print(f"E[y^2] = {vary}")
    print(f"Sigma = {sigy}")
    print(f"E[y^4] = {kury}")
    print(kury / (vary**2))

    # scatterplot()
    # print(kur_z)
    # print()
    # print(kur_z - Decimal(3))

    # print(f"\nExp x: {exp_x}\nSigma x: {sigma_x/Decimal(2)}\nVar x: {var_x}\nKur x: {kur_x}")
    # print(f"\nExp z: {exp_z}\nSigma z: {sigma_z}\nVar z: {var_z}\nKur z: {kur_z}")

    # num = (Decimal(0.68)**Decimal(2)) * Decimal(96)
    # den = (Decimal(10)**Decimal(-14))**Decimal(2)
    #
    # n = (num/den)
    # print(f"Required number of samples: {round(n)}")

    # num_samples = 1000000
    # t = [random.randint(0, 1) for _ in range(num_samples)]
    # x = sample_vector(t)
    # barchart(x)
    # # t = [0 for _ in range(num_samples)]
    # for i in range(-10, 10):
    #     print(f"{i}: {pr(i, 2 * 1.01)}")
    #
    #
    # print("T0")
    # for entry in Hawk256_T0_unscaled:
    #     print(f"{entry:.78f}")
    # print("T1")
    # for entry in Hawk256_T1_unscaled:
    #     print(f"{entry:.78f}")
    # print()
    # print("T0")
    # for entry in Hawk512_T0_unscaled:
    #     print(f"{entry:.78f}")
    # print("T1")
    # for entry in Hawk512_T1_unscaled:
    #     print(f"{entry:.78f}")
    # print()
    # print("T0")
    # for entry in Hawk1024_T0_unscaled:
    #     print(f"{entry:.78f}")
    # print("T1")
    # for entry in Hawk1024_T1_unscaled:
    #     print(f"{entry:.78f}")
