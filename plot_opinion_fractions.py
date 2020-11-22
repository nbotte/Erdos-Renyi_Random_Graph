import numpy as np
import matplotlib.pyplot as plt

fractions_av = np.loadtxt('Fraction_of_opinions_1_50_50_no_stubb_paper8_active_1_one_node_fully_connected_av100_long.txt')
#fractions_75_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_75_bern_050_av.txt')
#fractions_50_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_50_bern_050_av.txt')
#fractions_25_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_25_bern_050_av.txt')
#fractions_1_av = np.loadtxt('Fraction_of_opinions_1_70_30_all_stubb_100_bern_5_av.txt')

t = np.zeros(10000)
y = np.zeros(10000)

for i in range(len(t)):
    t[i] = i
    y[i] = 0.5

#plt.plot(t, fractions_0_av[:,0], label='Stubborness = 0')
#plt.plot(t[::50], fractions_av[:,0][::50], label='Opinion 0')
plt.plot(t[::50], fractions_av[:,1][::50], label='Opinion 1')
plt.plot(t[::50], y[::50])
#plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
plt.xlabel("Timesteps t")
plt.ylabel("Opinion fraction")
#plt.xlim(0, 50)
plt.ylim(0.45, 0.55)
plt.legend(loc='best')
plt.title('Opinion fraction vs time (probabilistic model, one node), 50/50, p = 1.\nFully connected network, averaged over 100 networks')
plt.savefig('Fraction_of_opinions_1_50_50_no_stubb_paper8_active_1_one_node_fully_av100_long_random_init.png')
plt.show()
