import numpy as np


def random_orthonormal_matrix(n):
    # Generate a random n x n matrix
    A = np.random.randn(n, n)
    # Perform QR decomposition
    Q, R = np.linalg.qr(A)
    return Q


def random_unit_vector(n):
    # Generate a random vector of size n
    v = np.random.randn(n)
    # Normalize the vector to have unit length
    v /= np.linalg.norm(v)
    return v


n = 16

C = random_orthonormal_matrix(n)
# u = C[:, 14]
# print(C)

avg = 0
samples = 1000
for i in range(samples):
    u = random_unit_vector(n)
    # print(np.linalg.norm())
    if np.sum(u.dot(C) ** 4) >= 1:
        print("Jasså")
    avg += np.sum(u.dot(C) ** 4)

avg /= samples
print()
print(avg)
