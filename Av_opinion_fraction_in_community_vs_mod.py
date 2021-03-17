import numpy as np
import matplotlib.pyplot as plt

highModPR = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=0-1_PR_01-0001_10x100_T=0.txt")
mediumModPR = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=0-1_PR_007-0004_10x100_T=0.txt")
lowModPR = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=0-1_PR_003-0008_10x100_T=0.txt")

highModPRrand = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=44-56_PR_01-0001_10x100_T=0_random=50-50.txt")
mediumModPRrand = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=44-56_PR_007-0004_10x100_T=0_random=50-50.txt")
lowModPRrand = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=44-56_PR_003-0008_10x100_T=0_random=50-50.txt")

highModREC = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=0-1_REC_01-0001_10x100_T=0.txt")
mediumModREC = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=0-1_REC_007-0004_10x100_T=0.txt")
lowModREC = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=0-1_REC_003-0008_10x100_T=0.txt")

highModRECrand = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=44-56_REC_01-0001_10x100_T=0_random=50-50_test.txt")
mediumModRECrand = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=44-56_REC_007-0004_10x100_T=0_random=50-50_test.txt")
lowModRECrand = np.loadtxt("Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=44-56_REC_003-0008_10x100_T=0_random=50-50_test.txt")

fracBegin = np.zeros(3)
fracEndPR = np.zeros(3)
fracEndPRrand = np.zeros(3)
fracEndPRsameComm = np.zeros(3)
fracEndREC = np.zeros(3)
fracEndRECsameComm = np.zeros(3)

fracBegin[0] = np.mean([np.mean(highModPR[:,0]), np.mean(highModREC[:,0]), np.mean(highModPRrand[:,0]), np.mean(highModRECrand[:,0])])
fracBegin[1] = np.mean([np.mean(mediumModPR[:,0]), np.mean(mediumModREC[:,0]), np.mean(mediumModPRrand[:,0]), np.mean(mediumModRECrand[:,0])])
fracBegin[2] = np.mean([np.mean(lowModPR[:,0]), np.mean(lowModREC[:,0]), np.mean(lowModPRrand[:,0]), np.mean(lowModRECrand[:,0])])

fracEndPR[0] = np.mean(highModPR[:,1])
fracEndPR[1] = np.mean(mediumModPR[:,1])
fracEndPR[2] = np.mean(lowModPR[:,1])

fracEndPRrand[0] = np.mean(highModPRrand[:,1])
fracEndPRrand[1] = np.mean(mediumModPRrand[:,1])
fracEndPRrand[2] = np.mean(lowModPRrand[:,1])

fracEndREC[0] = np.mean(highModREC[:,1])
fracEndREC[1] = np.mean(mediumModREC[:,1])
fracEndREC[2] = np.mean(lowModREC[:,1])

fracEndRECsameComm[0] = np.mean(highModRECrand[:,1])
fracEndRECsameComm[1] = np.mean(mediumModRECrand[:,1])
fracEndRECsameComm[2] = np.mean(lowModRECrand[:,1])

mod = np.zeros(3)
mod[0] = 0.9
mod[1] = 0.6
mod[2] = 0.25

y = [0.6, 0.6, 0.6]

randPR = [np.mean(fracEndPRrand), np.mean(fracEndPRrand), np.mean(fracEndPRrand)]


plt.plot(mod, randPR, 'g--')
#plt.plot(mod, y, 'r--')
plt.plot(mod, fracBegin, 'ro', label = r"Average fraction of opinion 0 at beginning, t = 0")
#plt.plot(mod, fracEndPR, 'bo', label = r"Average fraction of opinion 0 at end, t = 500; PR")
plt.plot(mod, fracEndPRrand, 'go', label = r"Average fraction of opinion 0 at end, t = 500; PR; 50/50")
#plt.plot(mod, fracEndREC, 'yo', label = r"Average fraction of opinion 0 at end, t = 500; REC")
plt.plot(mod, fracEndRECsameComm, 'co', label = r"Average fraction of opinion 0 at end, t = 500; REC; 50/50")
plt.xlabel("Modularity")
plt.ylabel("Average fraction of opinion 0 inside the communities")
#plt.ylim(0.5, 1.3)
plt.legend(loc='upper right')
plt.title("Opinion fraction inside communities before and after opinion evolution \n SBM (10 x 100) \n 0.5 communities with all nodes having opinion 0; rest 0/1; total 50/50")
plt.savefig("Fraction_of_opinion_0_inside_the_communities_before_and_after_the_opinion_evolution_vs_mod_SBM_10x100_3_diff_mod_comm=05_50-50.png")
plt.show()
