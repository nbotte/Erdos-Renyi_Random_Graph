import numpy as np
import matplotlib.pyplot as plt

clus_coeff = np.loadtxt('Clustering_coefficient_WS_vs_beta.txt')

plt.plot(clus_coeff[:,0], clus_coeff[:,1])
plt.show()
