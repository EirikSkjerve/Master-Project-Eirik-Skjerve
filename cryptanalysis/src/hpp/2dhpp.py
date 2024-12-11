import matplotlib.pyplot as plt
import numpy as np


def random_matrix_w_det_one():
    while True:
        # Generate random integers for the matrix entries
        matrix = np.random.randint(-2, 2, (2, 2))

        # Calculate the determinant
        det = np.linalg.det(matrix)

        # Check if the determinant is close to 1 (to account for floating-point precision)
        if round(det) == 1:
            return matrix


if __name__ == "__main__":
    # B = random_matrix_with_determinant_one()

    # B = np.array(([1, 1], [0, 1]))
    # B = random_matrix_w_det_one()
    B = np.array(([1,0],[0,1]))
    print(f"B: {B}")

    num_samples = 1000000
    mean = 0
    std_dev = 2

    zsn = []
    zsu = []
    for n in range(num_samples):

        zsn.append(np.random.normal(mean, std_dev, 2).round())
        zsu.append(np.random.uniform(-5, 5, 2).round())
        # zsn.append(np.random.normal(mean, std_dev, 2))
        # zsu.append(np.random.uniform(-5, 5, 2))

    wsn = [B@zsn[i] for i in range(num_samples)]
    wsu = [B@zsu[i] for i in range(num_samples)]

    xsn = [wsn[i][0] for i in range(num_samples)]
    ysn = [wsn[i][1] for i in range(num_samples)]

    xsu = [wsu[i][0] for i in range(num_samples)]
    ysu = [wsu[i][1] for i in range(num_samples)]

    # plt.grid(True)
    plt.scatter(xsn, ysn, color="green", alpha=0.1)
    plt.scatter(xsu, ysu, alpha=0.1)

    # # Generate points for the unit circle
    # theta = np.linspace(0, 2 * np.pi, 100)
    # x_circle = np.cos(theta)
    # y_circle = np.sin(theta)
    #
    # # Plot the unit circle on top
    # plt.plot(x_circle, y_circle, color='black', linewidth=2, label='Unit Circle')
    # plt.axhline(0, color='grey', lw=0.5)  # Add x-axis
    # plt.axvline(0, color='grey', lw=0.5)  # Add y-axis

    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
