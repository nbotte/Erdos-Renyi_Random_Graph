import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    return a * x + b

echo_chambers1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_10x100.txt")
echo_chambers2 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_20x50.txt")
echo_chambers3 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_50x20.txt")
echo_chambers4 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_100x10.txt")

echo = np.zeros(4)

echo[0] = echo_chambers1[:,2][0]
echo[1] = echo_chambers2[:,2][0]
echo[2] = echo_chambers3[:,2][0]
echo[3] = echo_chambers4[:,2][0]

sizes = np.zeros(4)

sizes[0] = 100
sizes[1] = 50
sizes[2] = 20
sizes[3] = 10

popt, pcov = curve_fit(func, sizes, np.log(echo))
print(popt)

plt.plot(sizes, np.log(echo), 'o', label = 'data points')
plt.plot(sizes, func(sizes, *popt), label="linear fit to log(y) = a*x + b")
plt.xlabel("Community sizes")
plt.ylabel("Fraction of nodes with all neigbors having the same opinion 0")
plt.legend(loc='best')
plt.title('Size of echo chamber versus community size, 50/50\n N = 1000, p_cl = 0.1, p_add = 0.001\nStochastic block model, 10 x 10 averaged')
#plt.savefig("echo_chamber_vs_community_size_log.png")
plt.show()
