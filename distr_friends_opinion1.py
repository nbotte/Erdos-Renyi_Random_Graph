import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_0_REC.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_0_PR.txt")

# probably need to change size of t and y back to 10!
t = np.zeros(7)
y = np.zeros(7)

for i in range(7):
    t[i] = i/7
    y[i] = 1

print(len(t))

plt.plot(t, hist[:,2], linestyle='solid', label="REC")
plt.plot(t, hist1[:,2], linestyle='solid', label="PR")
plt.plot(t, y)


plt.xlabel("Fraction of friends with opinion 1")
plt.ylabel("Normalized average distribution")
plt.ylim(0, 20)

plt.legend(loc='best')

plt.title('Normalized average distribution of friends with the same opinion 1\nbeta = 0, K = 6, p_act = 0.1, N = 1000\nWatts-Strogatz network, 10 x 10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_WS_0_REC_PR.png')
plt.show()
