import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_1-001_REC.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_1-001_PR.txt")

# probably need to change size of t and y back to 10!
t = np.zeros(10)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

print(len(t))

plt.plot(t, hist[:,2], linestyle='solid', label="REC")
plt.plot(t, hist1[:,2], linestyle='solid', label="PR")
plt.plot(t, y)


plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
plt.ylim(0, 20)

plt.legend(loc='best')

plt.title('Normalized average distribution of friends with the same opinion 1\np_cl = 0.1, p_add = 0.001\nStochastic block model (100 x 10), 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_SBM_1-001_100x10_REC_PR.png')
plt.show()
