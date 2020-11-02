import numpy as np
import matplotlib.pyplot as plt

#fractions_0_av = np.loadtxt('Fraction_of_opinions_1_70_30_no_stubb_bern_5_av.txt')
fractions_75_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_75_bern_050_av.txt')
fractions_50_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_50_bern_050_av.txt')
fractions_25_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_25_bern_050_av.txt')
#fractions_1_av = np.loadtxt('Fraction_of_opinions_1_70_30_all_stubb_100_bern_5_av.txt')

t = np.zeros(300)
y = np.zeros(300)

for i in range(len(t)):
    t[i] = i
    y[i] = 0.5

#plt.plot(t, fractions_0_av[:,0], label='Stubborness = 0')
plt.plot(t, fractions_75_av[:,0], label='Stubborness = 0.75')
plt.plot(t, fractions_50_av[:,0], label='Stubborness = 0.50')
plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
plt.xlabel("Timesteps t")
plt.ylabel("Opinion fraction of minority opinion")
plt.xlim(0, 50)
#plt.ylim(-0.01, 0.300)
plt.legend(loc='best')
plt.title('Minority opinion vs time (av. over 100 network), 80/20, p=0.1 \n Bernouilli probability for node to be active is 0.5\n 50% of nodes are stubborn, stubborness goes from 0.25 to 0.75')
plt.savefig('Fraction_of_opinions_1_80_20_50_stubb_bern_5_av.png')
plt.show()
