import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("heat_output.csv", delimiter=",")
plt.imshow(data, cmap='hot', interpolation='nearest')
plt.colorbar(label="Temperature")
plt.title("Heat Diffusion Simulation")
plt.show()
