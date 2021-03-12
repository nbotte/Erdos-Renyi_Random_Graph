import numpy as np
import matplotlib.pyplot as plt

fractions_rand = np.loadtxt('Fraction_of_opinions_SBM_active_01_av_good_init_commOp0=04_other=50-50_REC_01-0001_10x100_T=0_random=70-30.txt')
fractions_comm02 = np.loadtxt('Fraction_of_opinions_SBM_active_01_av_good_init_commOp0=04_other=50-50_REC_01-0001_10x100_T=0.txt')
fractions_sameComm02 = np.loadtxt('Fraction_of_opinions_SBM_active_01_av_good_init_commOp0=04_other=50-50_REC_01-0001_10x100_T=0_sameComm.txt')
#fractions_PR = np.loadtxt('Fraction_of_opinions_Clustered_Cluster5_01-0001_50_50_no_stubb_paper8_active_01_av_good_init_PR.txt')
#fractions_50_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_50_bern_050_av.txt')
#fractions_25_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_25_bern_050_av.txt')
#fractions_1_av = np.loadtxt('Fraction_of_opinions_1_70_30_all_stubb_100_bern_5_av.txt')

t = np.zeros(500)
y = np.zeros(500)

for i in range(len(t)):
    t[i] = i
    y[i] = 0.7

#plt.plot(t, fractions_0_av[:,0], label='Stubborness = 0')
#plt.plot(t[::10], fractions_av[:,0][::10], label='Opinion 0')
plt.plot(t[::10], fractions_rand[:,0][::10], label='70/30, REC')
plt.plot(t[::10], fractions_comm02[:,0][::10], label='0.4 community opinion 0, other 50/50; REC')
plt.plot(t[::10], fractions_sameComm02[:,0][::10], label='0.4 community opinion 0 (same comm.), other 50/50; REC')
plt.plot(t[::10], y[::10])
#plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
plt.xlabel("Timesteps t")
#plt.ylabel("Opinion fraction")
plt.ylabel("Opinion 0 fracion")
#plt.xlim(0, 50)
plt.ylim(0.5, 1.2)
plt.legend(loc='upper right')
plt.title('Opinion fraction vs time, 70/30\n' r'$10 x 100; p_{cl} = 0.1; p_{add} = 0.001$' '\nSBM, 10 x 10 averaged')
#plt.title('Standard deviation for each timestep\nStochastic block model, 10 x 10 averaged')
plt.savefig('Fraction_of_opinions_SBM_active_01_av_good_init_commOp0=04_other=50-50_REC_01-0001_10x100_T=0_random=70-30_sameComm.png')
#plt.savefig('Standard_Deviation_Clustered_01-001_50_50_no_stubb_one_node_active_1_av_good_init.png')

plt.show()

"""fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
st = fig.suptitle("Stochastic block model (100 x 10), 10 x 10 averaged")
#axes[0].plot(t, y)
#axes[0].plot(t[::10], fractions_av0[:,0][::10], label='Opinion 0')
axes[0].plot(t[::10], fractions_REC[:,1][::10], label='Opinion 1, REC')
axes[0].legend(loc='best')
axes[0].set_xlabel("Timesteps t")
axes[0].set_ylabel("Opinion fraction inside a community, REC")
axes[0].set_ylim(0.39, 0.55)
axes[0].title.set_text("Opinion fraction inside a community vs time, REC, 50/50\n"r"$p_{cl}$ = 0.1, $p_{add}$ = 0.001")

axes[1].plot(t[::10], fractions_PR[:,1][::10], label='Opinion 1, PR')
#axes[1].plot(t[::10], fractions_av3[:,1][::10], label='Standard deviation')
#axes[1].plot(t, y)
axes[1].legend(loc='best')
axes[1].set_xlabel("Timesteps t")
axes[1].set_ylabel("Opinion fraction inside a community, PR")
axes[1].set_ylim(0.40, 0.45)
axes[1].title.set_text("Opinion fraction inside a community vs time, PR, 50/50\n"r"$p_{cl}$ = 0.1, $p_{add}$ = 0.001")
# shift subplots down:
st.set_y(0.03)
fig.subplots_adjust(bottom=0.1)

fig.tight_layout()
fig.savefig('Fraction_of_opinions_Clustered_Cluster5_01-0001_50_50_no_stubb_paper8_active_01_av_good_init_100x10_REC_PR.png')
plt.show()"""
