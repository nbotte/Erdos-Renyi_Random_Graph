import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a*np.exp(b*x) + c

echo_chambers1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_10x100.txt")
echo_chambers2 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0005_res=0_10x100.txt")
echo_chambers3 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=0_10x100.txt")
echo_chambers4 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-005_res=0_10x100.txt")
echo_chambers5 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-01_res=0_10x100.txt")

echo_chambers11 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_10x100.txt")
echo_chambers21 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_20x50.txt")
echo_chambers31 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_50x20.txt")
echo_chambers41 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_100x10.txt")
echo_chambers51 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_200x5.txt")

echo = np.zeros(5)

echo[0] = echo_chambers1[:,2][0]
echo[1] = echo_chambers2[:,2][0]
echo[2] = echo_chambers3[:,2][0]
echo[3] = echo_chambers4[:,2][0]
echo[4] = echo_chambers5[:,2][0]

echo1 = np.zeros(5)

echo1[0] = echo_chambers11[:,2][0]
echo1[1] = echo_chambers21[:,2][0]
echo1[2] = echo_chambers31[:,2][0]
echo1[3] = echo_chambers41[:,2][0]
echo1[4] = echo_chambers51[:,2][0]

mod = np.zeros(5)

mod[0] = 0.910383
mod[1] = 0.655756
mod[2] = 0.476086
mod[3] = 0.098409
mod[4] = 0.00900886

mod1 = np.zeros(5)

mod1[0] = 0.910383
mod1[1] = 0.830371
mod1[2] = 0.651065
mod1[3] = 0.469297
mod1[4] = 0.286703

popt, pcov = curve_fit(func, mod, echo)
print(popt)
popt1, pcov1 = curve_fit(func, mod1, echo1)
print(popt1)
xf = np.linspace(0,0.92,50)

plt.plot(mod, echo, 'o', label = r'Changing $p_{add}$')
plt.plot(xf, func(xf, *popt))
plt.plot(mod1, echo1, 'o', label = r'Changing clustersizes')
plt.plot(xf, func(xf, *popt1))
plt.xlabel("Modularity")
plt.ylabel("Fraction of nodes with all neigbors having the same opinion 0")
plt.legend(loc='best')
plt.title('Size of echo chamber versus modularity, 50/50\n N = 1000 \nStochastic block model, 10 x 10 averaged')
plt.savefig("echo_chamber_vs_modularity_edge_prob_clustersizes.png")
plt.show()
