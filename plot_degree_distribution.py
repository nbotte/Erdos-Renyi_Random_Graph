import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a*x**(-b) + c

def log_norm(x, mu, sigma):
    return 1/(np.sqrt(2*np.pi)*sigma*x)*np.exp(-((np.log(x)-
   mu)**2)/(2*sigma**2))

deg_distr = np.loadtxt('Degree_distribution_real_network_lastfm_fin.txt')

N = 10
p = 0.01

k = np.arange(0, 20, 1)

y = deg_distr[:,1]/np.trapz(deg_distr[:,0], deg_distr[:,1])

for i in range(len(y)):
    print(y[i])

popt, pcov = curve_fit(log_norm, deg_distr[:,0], y)

x = np.linspace(1, 70, 1000)

#plt.plot(k, binom.pmf(k, N, p), '.', markersize=5, label = 'Binomial distribution \nN = 999, p = 0.01')
plt.plot(deg_distr[:,0], np.abs(y), 'o', markersize='5')
plt.plot(x, log_norm(x, *popt), 'r--')
plt.xlabel("Degree k")
plt.ylabel("Degree distribution for lastfm-fin network")
#plt.xlim(0, 1060)
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(loc='upper right')
plt.title('Degree distribution for lastfm-fin network \n N = 8003')
plt.savefig('lastfm_fin_degree_distribution.png')
plt.show()
