''' NOTE: START MAKING NICE PLOTS THAT YOU CAN INSERT IN THESIS'''

import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_007-002-0007_powerlaw_T=0.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_003-0008_10x100_T=0.txt")
hist2 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0008_100x10_T=0.txt")

# probably need to change size of t and y back to 10!
t = np.zeros(10)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

plt.plot(t, hist[:,2], linestyle='solid', label=r"SBM; powerlaw; $p_{cl<100} = 0.07, p_{cl\geqslant 100} = 0.02, p_{add} = 0.007$")
plt.plot(t, hist1[:,2], linestyle='solid', label=r"SBM; $10x100; p_{cl} = 0.03, p_{add} = 0.008$")
plt.plot(t, hist2[:,2], linestyle='solid', label=r"SBM; $100x10; p_{cl} = 0.3, p_{add} = 0.008$")
plt.plot(t, y)


plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
plt.ylim(0, 30)

plt.legend(loc='upper center')

plt.title('Normalized average distribution of friends with the same opinion 1 \n' r'N = 1000; PR method; T = 0; Deg $\sim 10$; Mod $\sim 0.2$' '\n 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_SBM_powerlaw_vs_regular_PR_lowMod.png')
plt.show()
