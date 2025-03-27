import numpy as np
import matplotlib.pyplot as plt


# Function to project a point onto a line defined by a point and a direction vector
def project(point, v):
    """Projects a point onto a line defined by point p0 and direction vector v."""
    v = v / np.linalg.norm(v)  # Normalize direction vector
    proj = np.dot(point, v) * v  # Projection formula
    return proj


num_vectors = 1000
n = 2

np.random.seed(1337)
vectors = np.random.uniform(-1, 1, (num_vectors, n)).round()
# vectors = np.random.normal(0, 1, (num_vectors, n)).round()

v = np.array([[3, 0.5], [0.5, 2]])

l = np.linalg.cholesky(np.linalg.inv(v.T @ v))

res = vectors @ v
res_cube = (vectors @ v) @ l
c = v @ l


# Step to plot in 2D space
# Extract the first two dimensions for plotting
x = res_cube[:, 0]  # First component
y = res_cube[:, 1]  # Second component

# Extract the first two columns of B for the arrows
v_arrows_x = c[:, 0]  # First column
v_arrows_y = c[:, 1]  # Second column

projections = np.array([project(p, v_arrows_x) for p in res_cube])

# Create the scatter plot
plt.figure(figsize=(8, 8))
plt.scatter(
    x, y, alpha=0.6, color="black", s=3, label="transformed signature sample"
)  # Adjust alpha for transparency
# Use quiver to plot the arrows from origin (0, 0)
# plt.quiver(0, 0, v_arrows_x[0], v_arrows_y[0], angles='xy', scale_units='xy', scale=1, color='blue', label='C Column 1')
# plt.quiver(0, 0, v_arrows_x[1], v_arrows_y[1], angles='xy', scale_units='xy', scale=1, color='red', label='C Column 2')
# plt.title('2D Scatter Plot of Vectors Transformed by Matrix V')
# plt.xlabel('')
# plt.ylabel('')

# Plot the line using its direction vector
t = np.linspace(-10, 10, 100)
line_points = np.array([t_i * v_arrows_x for t_i in t])  # Parametric equation

plt.plot(line_points[:, 0], line_points[:, 1], "k-", linewidth=2, label="Line")

# Draw arrows from each point to its projection using annotate()
for p, proj in zip(res_cube, projections):
    plt.annotate(
        "",
        xy=(proj[0], proj[1]),  # End at projection
        xytext=(p[0], p[1]),  # Start at original point
        arrowprops=dict(arrowstyle="->", color="gray", lw=1.5),
    )

# Plot projections
plt.scatter(projections[:, 0], projections[:, 1], color="red", s=5, label="Projections")

plt.xlim(-5, 5)  # Adjust x-axis limits as needed
plt.ylim(-5, 5)  # Adjust y-axis limits as needed
plt.grid()
plt.axhline(0, color="black", linewidth=0.5, ls="--")
plt.axvline(0, color="black", linewidth=0.5, ls="--")
plt.legend()
plt.show()
