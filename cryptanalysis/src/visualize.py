import numpy as np
import matplotlib. pyplot as plt

num_vectors = 5000
n = 2

vectors = np.random.uniform(-1, 1, (num_vectors, n))
vectors = np.random.normal(0, 1, (num_vectors, n))

v = np.array([[1, 0],
             [0, 1]])
# theta = np.radians(57)
# v = np.array([[2*np.cos(theta), 0],
#               [3.5*np.sin(theta), 4*np.cos(theta)]])
res = vectors @ v

# Step to plot in 2D space
# Extract the first two dimensions for plotting
x = res[:, 0].round()  # First component
y = res[:, 1].round()  # Second component

# Extract the first two columns of B for the arrows
v_arrows_x = v[:, 0]  # First column
v_arrows_y = v[:, 1]  # Second column


# Create the scatter plot
plt.figure(figsize=(8, 8))
plt.scatter(x, y, alpha=0.1, color='blue')  # Adjust alpha for transparency
# Use quiver to plot the arrows from origin (0, 0)
plt.quiver(0, 0, v_arrows_x[0], v_arrows_y[0], angles='xy', scale_units='xy', scale=1, color='green', label='V Column 1')
plt.quiver(0, 0, v_arrows_x[1], v_arrows_y[1], angles='xy', scale_units='xy', scale=1, color='purple', label='V Column 2')
plt.title('2D Scatter Plot of Vectors Transformed by Matrix V')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.xlim(-15, 15)  # Adjust x-axis limits as needed
plt.ylim(-15, 15)  # Adjust y-axis limits as needed
plt.grid()
plt.axhline(0, color='black',linewidth=0.5, ls='--')
plt.axvline(0, color='black',linewidth=0.5, ls='--')
plt.legend()
plt.show()
