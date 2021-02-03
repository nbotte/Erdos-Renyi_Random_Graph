import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a*np.exp(b*x) + c

echo_chambers1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_10x100.txt")
echo_chambers2 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=0_10x100.txt")
echo_chambers3 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_005-001_res=0_10x100.txt")
echo_chambers4 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_001-001_res=0_10x100.txt")
echo_chambers5 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_001-0001_res=0_10x100.txt")

echo = np.zeros(5)

echo[0] = echo_chambers1[:,2][0]
echo[1] = echo_chambers2[:,2][0]
echo[2] = echo_chambers3[:,2][0]
echo[3] = echo_chambers4[:,2][0]
echo[4] = echo_chambers5[:,2][0]

mod = np.zeros(5)

mod[0] = 0.910383
mod[1] = 0.476086
mod[2] = 0.291744
mod[3] = 0.00792583
mod[4] = 0.46936

popt, pcov = curve_fit(func, mod, echo)
print(popt)
xf = np.linspace(0,0.95,50)

plt.plot(mod, echo, 'o', label = 'data points')
plt.plot(xf, func(xf, *popt), label=r"exponential fit $y = ae^{bx} + c$")
plt.xlabel("Modularity")
plt.ylabel("Fraction of nodes with all neigbors having the same opinion 0")
plt.legend(loc='best')
plt.title('Size of echo chamber versus modularity, 50/50\n N = 1000 \nStochastic block model, 10 x 10 averaged')
plt.savefig("echo_chamber_vs_modularity_edge_prob.png")
plt.show()
