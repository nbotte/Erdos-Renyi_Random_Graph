import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_REC_beta-001_res-05.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_beta-001_res-05.txt")

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

plt.title('Normalized average distribution of friends with the same opinion 1\nK = 6, beta = 0.01, N = 1000, res = 0.5\nWatts-Strogatz model, 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_WS_REC_PR_beta-001_res-05.png')
plt.show()
