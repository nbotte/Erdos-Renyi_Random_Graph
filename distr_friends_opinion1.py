import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_Clustered_REC.txt")
#hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_001_REC_av.txt")

t = np.zeros(10)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

plt.plot(t, hist[:,2], label="REC")
#plt.plot(t, hist1[:,2], label="REC")
plt.plot(t, y)

plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
#plt.ylim(0, 15)

plt.legend(loc='best')

plt.title('Normalized average distribution of friends with the same opinion 1\np_cl = 0.5, p_add = 0.001, p_act = 0.1, N = 1000\nStochastic block model, 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_Clustered_REC_av.png')
plt.show()
