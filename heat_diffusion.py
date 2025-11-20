import numpy as np
import matplotlib.pyplot as plt

#Loads the final heat grid version for each file implementaiton
A1 = np.loadtxt("heat_output_seq.csv", delimiter=",") #Sequential version
A2 = np.loadtxt("heat_output_openmp.csv", delimiter=",") #Open MP version
A3 = np.loadtxt("heat_output_mpi.csv", delimiter=",") #MPI version

fig, axes = plt.subplots(1, 3, figsize=(15, 5)) #Creates the image with the three sided by side grids

#Loops through each dataset
for ax, data, title in zip(axes, [A1, A2, A3],
                           ['Sequential', 'OpenMP', 'MPI']):
    im = ax.imshow(data, cmap='hot', origin='lower') #Heat grid image
    ax.set_title(title) #Title for each subplot
    ax.axis('off') #Remove axis marks and labels for a easier to read image

fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label='Temperature (°C)') #Temperature scale for all plots
plt.suptitle('Heat Diffusion Comparison Across Implementations') 
plt.show() #show the image on the screen
fig.savefig('heat_diffusion_comparison.png') #save image as a png
