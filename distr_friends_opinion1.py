''' NOTE: START MAKING NICE PLOTS THAT YOU CAN INSERT IN THESIS'''

import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_commOp0=04_other=50-50_REC_01-0001_10x100_T=0.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_commOp0=04_other=50-50_REC_01-0001_10x100_T=0_random=70-30.txt")
hist2 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_commOp0=04_other=50-50_REC_003-0008_10x100_T=0.txt")
hist3 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_commOp0=04_other=50-50_REC_003-0008_10x100_T=0_random=70-30.txt")

# probably need to change size of t and y back to 10!
t = np.zeros(10)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

plt.plot(t, hist[:,2], linestyle='solid', label=r"$p_{cl} = 0.1, p_{add} = 0.001$; 0.4 comm. op 0, other 50/50; high mod")
plt.plot(t, hist1[:,2], linestyle='solid', label=r"$p_{cl} = 0.1, p_{add} = 0.001$; 70/30; high mod")
plt.plot(t, hist2[:,2], linestyle='solid', label=r"$p_{cl} = 0.03, p_{add} = 0.008$; 0.4 comm. op 0, other 50/50; low mod")
plt.plot(t, hist3[:,2], linestyle='solid', label=r"$p_{cl} = 0.03, p_{add} = 0.008$; 70/30; low mod")
plt.plot(t, y)


plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
#plt.ylim(0, 35)

plt.legend(loc='upper right')

plt.title('Normalized average distribution of friends with the same opinion 1 \n' r'SBM (10 x 100); N = 1000; REC method; T = 0; Deg $\sim 10$' '\n 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_SBM_REC_commOp0=04_other=50-50_10x100_T=0.png')
plt.show()
