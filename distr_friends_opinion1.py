import numpy as np
import matplotlib.pyplot as plt

hist = np.loadtxt("Normalized_hist_fraction_friends_opinion1_001_av.txt")

t = np.zeros(10)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

plt.plot(t, hist)
plt.plot(t, y)

plt.xlabel("Averaged fraction of friends with opinion 1")
plt.ylabel("Averaged normalized distribution")
plt.ylim(0, 2)

plt.title('Averaged normalized distribution of friends with the same opinion 1\np = 0.01, p_act = 0.1, N = 1000\nErdos-Renyi network, 10x10 averaged')
plt.savefig('Normalized_hist_fraction_friends_opinion1_001_av.png')
plt.show()
