''' NOTE: START MAKING NICE PLOTS THAT YOU CAN INSERT IN THESIS'''

import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_REC_01-0001_10x100_T=0.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM-WS_REC_10-001-0001_10x100_T=0.txt")
hist2 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_REC_10-001_T=0.txt")
hist3 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_REC_10-006_T=0.txt")

# probably need to change size of t and y back to 10!
t = np.zeros(10)
t1 = np.zeros(11)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

for i in range(11):
    t1[i] = i/11


plt.plot(t1, hist2[:,2], linestyle='solid', label=r"WS; $\beta = 0.01, K = 10; C \sim 0.6$")
plt.plot(t1, hist3[:,2], linestyle='solid', label=r"WS; $\beta = 0.06, K = 10; C \sim 0.55$")
plt.plot(t, hist1[:,2], linestyle='solid', label=r"SBM-WS; $10x100; \beta = 0.01, K = 10, p_{add} = 0.001; C \sim 0.55; mod \sim 0.9$")
plt.plot(t, hist[:,2], linestyle='solid', label=r"SBM; $10x100; p_{cl} = 0.1, p_{add} = 0.001; C \sim 0.08; mod \sim 0.9$")
plt.plot(t, y)


plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
plt.ylim(0, 35)

plt.legend(loc='upper right')

plt.title('Normalized average distribution of friends with the same opinion 1 \n' r'N = 1000; REC method; Deg $\sim 10$' '\n 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_SBM_WS_SBM-WS_REC.png')
plt.show()
