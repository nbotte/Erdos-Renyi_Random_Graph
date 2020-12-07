import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt

deg_distr = np.loadtxt('Degree_distribution_Watts-Strogatz_05_av.txt')

N = 999
p = 0.01

k = np.arange(0, 20, 1)

#plt.plot(k, binom.pmf(k, N, p), '.', markersize=5, label = 'Binomial distribution \nN = 999, p = 0.01')
plt.plot(deg_distr[:,0], deg_distr[:,1]/N, label = 'Average of 100 simulated \nWatts-Strogatz networks')
plt.xlabel("Degree k")
plt.ylabel("Average Watts-Strogatz degree distribution")
plt.xlim(0, 60)
plt.legend(loc='upper right')
plt.title('Average degree distribution for Watts-Strogatz network \n N = 1000, K = 20, ' r'$\beta$ = 0.5, averaged over 100 independent networks')
plt.savefig('WS_degree_distribution_05_av.png')
plt.show()
