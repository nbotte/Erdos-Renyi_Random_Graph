import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt

deg_distr = np.loadtxt('Degree_distribution_Erdos-Renyi_01_av.txt')

N = 999
p = 0.01

k = np.arange(0, 20, 1)

plt.plot(k, binom.pmf(k, N, p), '.', markersize=5, label = 'Binomial distribution \nN = 999, p = 0.01')
plt.plot(deg_distr[:,0], deg_distr[:,1]/N, label = 'Average of 100 simulated \nErdos-Renyi network')
plt.xlabel("Degree k")
plt.ylabel("Average Erdos-Renyi degree distribution")
plt.xlim(0, 20)
plt.legend(loc='upper right')
plt.title('Average degree distribution for Erdos-Renyi network \n N = 1000, p = 0.01, averaged over 100 independent networks')
plt.savefig('ER_degree_distribution_01_av.png')
plt.show()
