import numpy as np
import matplotlib.pyplot as plt

fractions_av0 = np.loadtxt('Fraction_of_opinions_Clustered_01-001_50_50_no_stubb_paper8_active_01_av_good_init.txt')
#fractions_av3 = np.loadtxt('Fraction_of_opinions_Clustered_cluster6_01-001_50_50_no_stubb_paper8_active_01_av_good_init.txt')
#fractions_50_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_50_bern_050_av.txt')
#fractions_25_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_25_bern_050_av.txt')
#fractions_1_av = np.loadtxt('Fraction_of_opinions_1_70_30_all_stubb_100_bern_5_av.txt')

t = np.zeros(500)
y = np.zeros(500)

for i in range(len(t)):
    t[i] = i
    y[i] = 0.5

#plt.plot(t, fractions_0_av[:,0], label='Stubborness = 0')
#plt.plot(t[::10], fractions_av[:,0][::10], label='Opinion 0')
"""plt.plot(t[::50], fractions_av[:,2][::50], label='Standard deviation')
#plt.plot(t[::50], y[::50])
#plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
plt.xlabel("Timesteps t")
#plt.ylabel("Opinion fraction")
plt.ylabel("Standard deviation")
#plt.xlim(0, 50)
#plt.ylim(0.45, 0.55)
plt.legend(loc='best')
#plt.title('Opinion fraction vs time (probabilistic model), 50/50\n p_cl = 0.1, p_add = 0.01\nStochastic block model, 10 x 10 averaged')
plt.title('Standard deviation for each timestep\nStochastic block model, 10 x 10 averaged')
#plt.savefig('Fraction_of_opinions_Clustered_01-001_50_50_no_stubb_one_node_active_1_av_good_init.png')
plt.savefig('Standard_Deviation_Clustered_01-001_50_50_no_stubb_one_node_active_1_av_good_init.png')

plt.show()"""

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
st = fig.suptitle("Stochastic block model, 10 x 10 averaged")
axes[0].plot(t, y)
axes[0].plot(t[::10], fractions_av0[:,0][::10], label='Opinion 0')
#axes[0].plot(t[::10], fractions_av0[:,1][::10], label='Opinion 1')
axes[0].legend(loc='best')
axes[0].set_xlabel("Timesteps t")
axes[0].set_ylabel("Opinion fraction")
axes[0].set_ylim(0.45, 0.55)
axes[0].title.set_text("Opinion fraction vs time, 50/50\n p_cl = 0.5, p_add = 0.001, p_act = 0.1")

axes[1].plot(t[::10], fractions_av0[:,2][::10], label='Standard deviation')
#axes[1].plot(t[::10], fractions_av3[:,1][::10], label='Standard deviation')
#axes[1].plot(t, y)
axes[1].legend(loc='best')
axes[1].set_xlabel("Timesteps t")
axes[1].set_ylabel("Standard deviation")
#axes[1].set_ylim(0.35, 0.55)
axes[1].title.set_text("Standard deviation for each timestep\n p_cl = 0.5, p_add = 0.001, p_act = 0.1")
# shift subplots down:
st.set_y(0.03)
fig.subplots_adjust(bottom=0.1)

fig.tight_layout()
fig.savefig('Fraction_of_opinions_Clustered_05-0001_50_50_no_stubb_paper8_active_01_av_good_init.png')
plt.show()
