import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt

deg_distr = np.loadtxt('Degree distribution Erdos-Renyi.txt')

N = 999
p = 0.1

k = np.arange(0, 200, 1)

plt.plot(k, binom.pmf(k, N, p), '.', markersize=5, label = 'Binomial distribution, N = 999, p = 0.1')
plt.plot(deg_distr[:,0], deg_distr[:,1]/N, label = 'Simulated Erdos-Renyi network')
plt.xlabel("Degree k")
plt.ylabel("Erdos-Renyi degree distribution")
plt.legend(loc='best')
plt.title('Degree distribution for Erdos-Renyi network with N = 1000, p = 0.1')
plt.savefig('ER_degree_distribution.png')
plt.show()
