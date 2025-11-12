import numpy as np
import matplotlib.pyplot as plt

A1 = np.loadtxt("heat_output_seq.csv", delimiter=",")
A2 = np.loadtxt("heat_output_openmp.csv", delimiter=",")
A3 = np.loadtxt("heat_output_mpi.csv", delimiter=",")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, data, title in zip(axes, [A1, A2, A3],
                           ['Sequential', 'OpenMP', 'MPI']):
    im = ax.imshow(data, cmap='hot', origin='lower')
    ax.set_title(title)
    ax.axis('off')

fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label='Temperature (°C)')
plt.suptitle('Heat Diffusion Comparison Across Implementations')
plt.show()
fig.savefig('heat_diffusion_comparison.png')