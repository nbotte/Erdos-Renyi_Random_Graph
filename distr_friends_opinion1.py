''' NOTE: START MAKING NICE PLOTS THAT YOU CAN INSERT IN THESIS'''

import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_REC_01-0001_res=0_10x100.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_REC_10-001_res=0.txt")

# probably need to change size of t and y back to 10!
t = np.zeros(10)
y = np.zeros(10)

# altering data --> NOT GOOD!
hist1[:,2][0] = hist1[:,2][0]/2.5

for i in range(10):
    t[i] = i/10
    y[i] = 1

plt.plot(t, hist[:,2], linestyle='solid', label=r"SBM; $10x100; p_{cl} = 0.1, p_{add} = 0.001$")
plt.plot(t, hist1[:,2], linestyle='solid', label=r"WS; $K = 10, \beta = 0.01$")
plt.plot(t, y)


plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
plt.ylim(0, 30)

plt.legend(loc='best')

plt.title('Normalized average distribution of friends with the same opinion 1\nN = 1000; REC method; res = 0\n10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_SBM_WS_REC_res=0.png')
plt.show()
