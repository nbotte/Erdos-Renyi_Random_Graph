import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

'''echoSBMH0 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=0.txt")
echoSBMH01 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=01.txt")
echoSBMH02 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=02.txt")
echoSBMH03 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=03.txt")
echoSBMH04 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=04.txt")
echoSBMH05 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=05.txt")
echoSBMH06 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=06.txt")
echoSBMH07 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=07.txt")
echoSBMH08 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=08.txt")
echoSBMH09 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=09.txt")
echoSBMH1 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=1.txt")

echoSBM0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=0_50x20.txt")
echoSBM01 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=01_50x20.txt")
echoSBM02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=02_50x20.txt")
echoSBM03 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=03_50x20.txt")
echoSBM04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=04_50x20.txt")
echoSBM05 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=05_50x20.txt")
echoSBM06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=06_50x20.txt")

echoSBML0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=0_100x10.txt")
echoSBML01 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=01_100x10.txt")
echoSBML02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=02_100x10.txt")
echoSBML03 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=03_100x10.txt")
echoSBML04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=04_100x10.txt")
echoSBML05 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=05_100x10.txt")
echoSBML06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=06_100x10.txt")

echoER0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=0.txt")
echoER02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=02.txt")
echoER04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=04.txt")
echoER06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=06.txt")'''

echoWS0 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=0.txt")
echoWS01 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=01.txt")
echoWS02 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=02.txt")
echoWS03 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=03.txt")
echoWS04 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=04.txt")
echoWS05 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=05.txt")
echoWS06 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=06.txt")
echoWS07 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=07.txt")
echoWS08 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=08.txt")
echoWS09 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=09.txt")
echoWS1 = np.loadtxt("Echo_chamber_WS_PR_6_001_T=1.txt")

def func(x, a, b, c):
    return a*np.exp(b*x) + c

def func_lin(x, a, b):
    return a*x + b

threshold = np.zeros(11)
for i in range(11):
    threshold[i] = i/10

echoWS = np.zeros(11)

echoWS[0] = np.mean(echoWS0[:,2])
echoWS[1] = np.mean(echoWS01[:,2])
echoWS[2] = np.mean(echoWS02[:,2])
echoWS[3] = np.mean(echoWS03[:,2])
echoWS[4] = echoWS04[:,2][0]
echoWS[5] = np.mean(echoWS05[:,2])
echoWS[6] = np.mean(echoWS06[:,2])
echoWS[7] = np.mean(echoWS07[:,2])
echoWS[8] = np.mean(echoWS08[:,2])
echoWS[9] = np.mean(echoWS09[:,2])
echoWS[10] = np.mean(echoWS1[:,2])

'''echoSBMH = np.zeros(11)

echoSBMH[0] = np.mean(echoSBMH0[:,2])
echoSBMH[1] = np.mean(echoSBMH01[:,2])
echoSBMH[2] = np.mean(echoSBMH02[:,2])
echoSBMH[3] = np.mean(echoSBMH03[:,2])
echoSBMH[4] = np.mean(echoSBMH04[:,2])
echoSBMH[5] = np.mean(echoSBMH05[:,2])
echoSBMH[6] = np.mean(echoSBMH06[:,2])
echoSBMH[7] = np.mean(echoSBMH07[:,2])
echoSBMH[8] = np.mean(echoSBMH08[:,2])
echoSBMH[9] = np.mean(echoSBMH09[:,2])
echoSBMH[10] = np.mean(echoSBMH1[:,2])'''



poptWS, pcovWS = curve_fit(func_lin, threshold[:6], echoWS[:6])
print(poptWS)
poptWS1, pcovWS1 = curve_fit(func, threshold[6:], echoWS[6:])
print(poptWS1)
xf = np.linspace(0,0.6,50)
xf1 = np.linspace(0.6,1,50)

plt.plot(threshold, echoWS, 'ro', label = r'WS, $K = 6; \beta = 0.01$')
plt.plot(xf, func_lin(xf, *poptWS), 'r--')
plt.plot(xf1, func(xf1, *poptWS1), 'r--')
#plt.plot(threshold, echoSBMH, 'bo', label = r'SBM, 10x100, $p_{cl} = 0.1; p_{add} = 0.001$, high mod')

plt.xlabel("Threshold to change opinion")
plt.ylabel("Fraction of nodes with all neigbors having the same opinion (average of 0 and 1)\n(echo chamber size)")
plt.legend(loc='best')
plt.title('Size of echo chamber versus threshold to change opinion, 50/50\n N = 1000, PR method \n10 x 10 averaged')
plt.savefig("echo_chamber_vs_threshold_SBM_ER_WS_PR.png")
plt.show()
