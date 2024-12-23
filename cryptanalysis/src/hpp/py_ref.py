'''

Reference HPP attack by Martin Feussner

'''

import numpy as np
import sympy as sp
import time
import random

# optional
np.set_printoptions(precision=2)


def generate_secret_matrix(N, entry_bound):
    matrix = []
    rank = 0

    while rank<N:
        matrix = np.random.randint(-entry_bound, entry_bound+1, size=(N, N))
        rank = np.linalg.matrix_rank(matrix)

    return matrix


def generate_random_unit_vector(n):
    random_vector = np.random.uniform(-1, 1, n)
    w = random_vector / np.linalg.norm(random_vector)
    return w


def approx_nabla_mom4(U, w):
    uw = np.dot(U, w)
    uw3 = uw ** 3
    uw3u = uw3[:, np.newaxis] * U
    g = 4 * np.sum(uw3u, axis=0) / U.shape[0]

    return g


def approx_mom4(U, w):
    uw = np.dot(U, w)
    uw4 = uw ** 4
    m = np.mean(uw4)
    return m


def HPP_gradient_descent(U, L_inverse, lr):
    N = U.shape[1]
    solutions = set()

    iterations = 0

    while len(solutions) < N:
        # Step 1
        w = generate_random_unit_vector(N)

        while True:
            iterations += 1
            # Step 2
            g = approx_nabla_mom4(U, w)

            # Step 3
            w_new = w - lr*g
            w_new = w_new / np.linalg.norm(w_new)

            # Step 4
            m = approx_mom4(U, w)
            m_new = approx_mom4(U, w_new)

            if m_new >= m:
                v = w @ L_inverse
                v = np.rint(v).astype(int)
                nv = -v

                if tuple(v) not in solutions and tuple(nv) not in solutions:
                    solutions.add(tuple(v))
                break
            else:
                w = w_new

    print(f"gradient descent used {iterations} iterations")
    return np.array(list(solutions), dtype=int)


def map_matching_rows(V_approximation, V):
    for i in range(V.shape[0]):
        for j in range(V.shape[0]):
            if np.array_equal(V_approximation[i,:], V[j,:]) or np.array_equal(-V_approximation[i,:], V[j,:]):
                print(f"{i+1} -> {j+1}")


if __name__ == '__main__':

    # Set a seed and some parameters
    np.random.seed(random.randint(0, 10000))
    N = 16
    entry_bound = 1
    dist_bound = 1
    num_samples = 5000 # preprocessing experiments to find a suitable number

    # Generate invertible secret matrix V
    V = generate_secret_matrix(N, entry_bound)

    # Getting the expectation of x^2
    x = sp.symbols('x')
    f = (x ** 2) / (2*dist_bound)
    Ex2 = float(sp.integrate(f, (x, -dist_bound, dist_bound)))

    # Collect samples P(V)
    Y = np.zeros((num_samples, N))
    for i in range(num_samples):
        x = np.random.uniform(-dist_bound, dist_bound, size=N)
        print(f"x: {x}")
        y = x @ V
        Y[i, :] = y

    # Compute approximation of Gram Matrix
    G_approximate = np.round((1/Ex2) * (Y.T@Y)/num_samples).astype(int)
    G_approximate_inverse = np.linalg.inv(G_approximate)
    # print(f"G appinv: {G_approximate_inverse}")
    L = np.linalg.cholesky(G_approximate_inverse)
    L_inverse = np.linalg.inv(L)

    # Samples of P(C) where C=VL, which is x@C
    U = Y@L / dist_bound
    learning_rate = 0.75 # delta -> descent parameter

    start_time = time.time()
    V_approximation = HPP_gradient_descent(U, L_inverse, learning_rate)
    end_time = time.time()

    print()
    print(f"Time taken for HPP_gradient_descent ... {end_time - start_time:.2f} seconds")
    print()
    print("Approximation of rows of +-V:")
    print(V_approximation)
    print()
    print("Actual V:")
    print(V)
    print()
    map_matching_rows(V_approximation, V)
