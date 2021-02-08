''' NOTE: START MAKING NICE PLOTS THAT YOU CAN INSERT IN THESIS'''

import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_resComm=05_10x100.txt")
#hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_beta-001_res-05.txt")

# probably need to change size of t and y back to 10!
t = np.zeros(10)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

plt.plot(t, hist[:,2], linestyle='solid', label="PR")
#plt.plot(t, hist1[:,2], linestyle='solid', label="PR")
plt.plot(t, y)


plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
plt.ylim(0, 20)

plt.legend(loc='best')

plt.title('Normalized average distribution of friends with the same opinion 1\nN = 1000, p_cl = 0.1, p_add = 0.001; res = 0.5\nStochastic block model (10 x 100), 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_SBM_PR_01-0001_resComm=05_10x100.png')
plt.show()
