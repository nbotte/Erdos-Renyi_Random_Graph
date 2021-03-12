import numpy as np
import matplotlib.pyplot as plt

opPRhigh01 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=50-50_PR_01-0001_10x100_T=0.txt")
opPRhigh02 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=02_other=50-50_PR_01-0001_10x100_T=0.txt")
opPRhigh03 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=03_other=50-50_PR_01-0001_10x100_T=0.txt")
opPRhigh04 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=04_other=50-50_PR_01-0001_10x100_T=0.txt")
opPRhigh05 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=50-50_PR_01-0001_10x100_T=0.txt")
opPRhigh06 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=06_other=50-50_PR_01-0001_10x100_T=0.txt")

opPRrandhigh01 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=50-50_PR_01-0001_10x100_T=0_random=55-45.txt")
opPRrandhigh02 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=02_other=50-50_PR_01-0001_10x100_T=0_random=60-40.txt")
opPRrandhigh03 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=03_other=50-50_PR_01-0001_10x100_T=0_random=65-35.txt")
opPRrandhigh04 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=04_other=50-50_PR_01-0001_10x100_T=0_random=70-30.txt")
opPRrandhigh05 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=50-50_PR_01-0001_10x100_T=0_random=75-25.txt")
opPRrandhigh06 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=06_other=50-50_PR_01-0001_10x100_T=0_random=80-20.txt")

opPRlow01 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=50-50_PR_003-0008_10x100_T=0.txt")
opPRlow02 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=02_other=50-50_PR_003-0008_10x100_T=0.txt")
opPRlow03 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=03_other=50-50_PR_003-0008_10x100_T=0.txt")
opPRlow04 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=04_other=50-50_PR_003-0008_10x100_T=0.txt")
opPRlow05 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=50-50_PR_003-0008_10x100_T=0.txt")
opPRlow06 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=06_other=50-50_PR_003-0008_10x100_T=0.txt")

opPRrandlow01 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=50-50_PR_003-0008_10x100_T=0_random=55-45.txt")
opPRrandlow02 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=02_other=50-50_PR_003-0008_10x100_T=0_random=60-40.txt")
opPRrandlow03 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=03_other=50-50_PR_003-0008_10x100_T=0_random=65-35.txt")
opPRrandlow04 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=04_other=50-50_PR_003-0008_10x100_T=0_random=70-30.txt")
opPRrandlow05 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=50-50_PR_003-0008_10x100_T=0_random=75-25.txt")
opPRrandlow06 = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=06_other=50-50_PR_003-0008_10x100_T=0_random=80-20.txt")


frac = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

opBegin = np.zeros(6)

opBegin[0] = np.mean([np.mean(opPRhigh01[:,0]), np.mean(opPRrandhigh01[:,0]), np.mean(opPRlow01[:,0]), np.mean(opPRrandlow01[:,0])])
opBegin[1] = np.mean([np.mean(opPRhigh02[:,0]), np.mean(opPRrandhigh02[:,0]), np.mean(opPRlow02[:,0]), np.mean(opPRrandlow02[:,0])])
opBegin[2] = np.mean([np.mean(opPRhigh03[:,0]), np.mean(opPRrandhigh03[:,0]), np.mean(opPRlow03[:,0]), np.mean(opPRrandlow03[:,0])])
opBegin[3] = np.mean([np.mean(opPRhigh04[:,0]), np.mean(opPRrandhigh04[:,0]), np.mean(opPRlow04[:,0]), np.mean(opPRrandlow04[:,0])])
opBegin[4] = np.mean([np.mean(opPRhigh05[:,0]), np.mean(opPRrandhigh05[:,0]), np.mean(opPRlow05[:,0]), np.mean(opPRrandlow05[:,0])])
opBegin[5] = np.mean([np.mean(opPRhigh06[:,0]), np.mean(opPRrandhigh06[:,0]), np.mean(opPRlow06[:,0]), np.mean(opPRrandlow06[:,0])])

opPRhigh = np.zeros(6)

opPRhigh[0] = np.mean(opPRhigh01[:,1])
opPRhigh[1] = np.mean(opPRhigh02[:,1])
opPRhigh[2] = np.mean(opPRhigh03[:,1])
opPRhigh[3] = np.mean(opPRhigh04[:,1])
opPRhigh[4] = np.mean(opPRhigh05[:,1])
opPRhigh[5] = np.mean(opPRhigh06[:,1])

opPRrandhigh = np.zeros(6)

opPRrandhigh[0] = np.mean(opPRrandhigh01[:,1])
opPRrandhigh[1] = np.mean(opPRrandhigh02[:,1])
opPRrandhigh[2] = np.mean(opPRrandhigh03[:,1])
opPRrandhigh[3] = np.mean(opPRrandhigh04[:,1])
opPRrandhigh[4] = np.mean(opPRrandhigh05[:,1])
opPRrandhigh[5] = np.mean(opPRrandhigh06[:,1])

opPRlow = np.zeros(6)

opPRlow[0] = np.mean(opPRlow01[:,1])
opPRlow[1] = np.mean(opPRlow02[:,1])
opPRlow[2] = np.mean(opPRlow03[:,1])
opPRlow[3] = np.mean(opPRlow04[:,1])
opPRlow[4] = np.mean(opPRlow05[:,1])
opPRlow[5] = np.mean(opPRlow06[:,1])

opPRrandlow = np.zeros(6)

opPRrandlow[0] = np.mean(opPRrandlow01[:,1])
opPRrandlow[1] = np.mean(opPRrandlow02[:,1])
opPRrandlow[2] = np.mean(opPRrandlow03[:,1])
opPRrandlow[3] = np.mean(opPRrandlow04[:,1])
opPRrandlow[4] = np.mean(opPRrandlow05[:,1])
opPRrandlow[5] = np.mean(opPRrandlow06[:,1])

plt.plot(frac, opBegin, 'k', label = 'Average fraction in beginning, t = 0')
plt.plot(frac, opPRrandhigh, '--', label = r'Evenly distributed; high mod ($p_{cl} = 0.1, p_{add} = 0.001$)')
plt.plot(frac, opPRhigh, '--', label = r'Communities with single opinion, others 50/50; high mod')
plt.plot(frac, opPRrandlow, '--', label = r'Evenly distributed; low mod ($p_{cl} = 0.03, p_{add} = 0.008$)')
plt.plot(frac, opPRlow, '--', label = r'Communities with single opinion, others 50/50; low mod')

plt.xlabel("Fraction of communities with a single opinion equal to 0")
plt.ylabel("Average fraction of opinion 0 inside the communities")
plt.ylim(0.5, 1.3)
plt.title("Average fraction of opinion 0 vs fraction communities with single opinion \n SBM (10 x 100), PR method")
plt.legend(loc='best')
plt.savefig("fraction_of_opinion_0_vs_community_fraction_SBM_10x100_PR_high-low_mod.png")

plt.show()
