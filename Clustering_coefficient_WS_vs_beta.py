import numpy as np
import matplotlib.pyplot as plt

clus_coeff = np.loadtxt('Clustering_coefficient_WS_vs_beta.txt')

def func(beta):
    return (1-beta)**3

plt.plot(clus_coeff[:,0], clus_coeff[:,1]*4*19/(3*18), label = r'$\frac{C(\beta)}{C(0)}$')
plt.plot(clus_coeff[:,0], clus_coeff[:,1], label = r'C($\beta$)')
plt.plot(clus_coeff[:,0], func(clus_coeff[:,0]), label = r'$(1-\beta)^3$')
plt.legend(loc='best')
plt.xlabel(r'Rewire probability $\beta$')
plt.ylabel('Clustering coefficient')
plt.title(r'Clustering coefficient C vs $\beta$')
plt.savefig('Clustering_coefficient_WS_vs_beta.png')
plt.show()
