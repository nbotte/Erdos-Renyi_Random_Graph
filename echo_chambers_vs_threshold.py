import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

'''echoSBMH0 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=0.txt")
echoSBMH01 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=01.txt")
echoSBMH02 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=02.txt")
echoSBMH03 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=03.txt")
echoSBMH04 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=04.txt")
echoSBMH05 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=05.txt")
echoSBMH06 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=06.txt")
echoSBMH07 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=07.txt")
echoSBMH08 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=08.txt")
echoSBMH09 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=09.txt")
echoSBMH1 = np.loadtxt("Echo_chamber_SBM_REC_007-00005_10x100_T=1.txt")

echoSBM0 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=0.txt")
echoSBM01 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=01.txt")
echoSBM02 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=02.txt")
echoSBM03 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=03.txt")
echoSBM04 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=04.txt")
echoSBM05 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=05.txt")
echoSBM06 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=06.txt")
echoSBM07 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=07.txt")
echoSBM08 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=08.txt")
echoSBM09 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=09.txt")
echoSBM1 = np.loadtxt("Echo_chamber_SBM_REC_025-00025_50x20_T=1.txt")

echoSBML0 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=0.txt")
echoSBML01 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=01.txt")
echoSBML02 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=02.txt")
echoSBML03 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=03.txt")
echoSBML04 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=04.txt")
echoSBML05 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=05.txt")
echoSBML06 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=06.txt")
echoSBML07 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=07.txt")
echoSBML08 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=08.txt")
echoSBML09 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=09.txt")
echoSBML1 = np.loadtxt("Echo_chamber_SBM_REC_006-0007_100x10_T=1.txt")'''

'''echoER0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=0.txt")
echoER02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=02.txt")
echoER04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=04.txt")
echoER06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=06.txt")'''

echoWS0 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=0.txt")
echoWS01 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=01.txt")
echoWS02 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=02.txt")
echoWS03 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=03.txt")
echoWS04 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=04.txt")
echoWS05 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=05.txt")
echoWS06 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=06.txt")
echoWS07 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=07.txt")
echoWS08 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=08.txt")
echoWS09 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=09.txt")
echoWS1 = np.loadtxt("Echo_chamber_WS_REC_6-001_T=1.txt")

def func(x, a, b, c):
    return a*np.exp(b*x) + c

def func_lin(x, a, b):
    return a*x + b

def func_trip(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

threshold = np.zeros(11)
y = np.zeros(11)
for i in range(11):
    threshold[i] = i/10
    y[i] = 1

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
echoSBMH[10] = np.mean(echoSBMH1[:,2])

echoSBM = np.zeros(11)

echoSBM[0] = np.mean(echoSBM0[:,2])
echoSBM[1] = echoSBM01[:,2][0]
echoSBM[2] = echoSBM02[:,2][-1]
echoSBM[3] = np.mean(echoSBM03[:,2])
echoSBM[4] = np.mean(echoSBM04[:,2])
echoSBM[5] = np.mean(echoSBM05[:,2])
echoSBM[6] = np.mean(echoSBM06[:,2])
echoSBM[7] = np.mean(echoSBM07[:,2])
echoSBM[8] = np.mean(echoSBM08[:,2])
echoSBM[9] = np.mean(echoSBM09[:,2])
echoSBM[10] = np.mean(echoSBM1[:,2])

echoSBML = np.zeros(11)

echoSBML[0] = echoSBML0[:,2][-1]
echoSBML[1] = echoSBML01[:,2][0]
echoSBML[2] = echoSBML02[:,2][-1]
echoSBML[3] = np.mean(echoSBML03[:,2])
echoSBML[4] = np.mean(echoSBML04[:,2])
echoSBML[5] = np.mean(echoSBML05[:,2])
echoSBML[6] = np.mean(echoSBML06[:,2])
echoSBML[7] = np.mean(echoSBML07[:,2])
echoSBML[8] = np.mean(echoSBML08[:,2])
echoSBML[9] = np.mean(echoSBML09[:,2])
echoSBML[10] = np.mean(echoSBML1[:,2])'''



poptWS, pcovWS = curve_fit(func_lin, threshold[:7], echoWS[:7])
print(poptWS)
poptWS1, pcovWS1 = curve_fit(func_trip, threshold[7:], echoWS[7:])
print(poptWS1)
xfWS = np.linspace(0,0.7,50)
xfWS1 = np.linspace(0.7,1,50)
WSerr = [5.4, 6.9, 6.7, 6.3, 6.3, 5.9, 5.8, 7.1, 3.1, 1.4, 0.5]

'''poptSBMH, pcovSBMH = curve_fit(func_lin, threshold[:5], echoSBMH[:5])
print(poptSBMH)
poptSBMH1, pcovSBMH1 = curve_fit(func_trip, threshold[5:], echoSBMH[5:])
print(poptSBMH1)
xfSBMH = np.linspace(0,0.5,50)
xfSBMH1 = np.linspace(0.5,1,50)
SBMHerr = [2.8, 3, 2.2, 2.7, 2.6, 2.6, 1.5, 0.9, 0.5, 0.4, 0.3]

poptSBM, pcovSBM = curve_fit(func_lin, threshold[:5], echoSBM[:5])
print(poptSBM)
poptSBM1, pcovSBM1 = curve_fit(func_trip, threshold[5:], echoSBM[5:])
print(poptSBM1)
xfSBM = np.linspace(0,0.5,50)
xfSBM1 = np.linspace(0.5,1,50)
SBMerr = [1.8, 1.7, 1.3, 1.2, 1.7, 1.7, 0.9, 0.6, 0.5, 0.4, 0.3]

poptSBML, pcovSBML = curve_fit(func_lin, threshold[:5], echoSBML[:5])
print(poptSBML)
poptSBML1, pcovSBML1 = curve_fit(func_trip, threshold[5:], echoSBML[5:])
print(poptSBML1)
xfSBML = np.linspace(0,0.5,50)
xfSBML1 = np.linspace(0.5,1,50)
SBMLerr = [0.9, 1, 1.3, 1.4, 1.6, 1, 0.7, 0.5, 0.4, 0.3, 0.3]'''

plt.plot(threshold, echoWS, 'ro', label = r'WS, $K = 6; \beta = 0.01$')
plt.plot(xfWS, func_lin(xfWS, *poptWS), 'r--')
plt.plot(xfWS1, func_trip(xfWS1, *poptWS1), 'r--')
#plt.errorbar(threshold, echoWS, WSerr, c='r', fmt='o', label = r'WS, $K = 6; \beta = 0.01$')

#plt.plot(threshold, echoSBMH, 'bo', label = r'SBM, 10x100, $p_{cl} = 0.07; p_{add} = 0.0005$, high mod')
'''plt.plot(xfSBMH, func_lin(xfSBMH, *poptSBMH), 'b--')
plt.plot(xfSBMH1, func_trip(xfSBMH1, *poptSBMH1), 'b--')
plt.errorbar(threshold, echoSBMH, SBMHerr, c='b', fmt='o', label = r'SBM, 10x100, $p_{cl} = 0.07; p_{add} = 0.0005$, high mod')

plt.plot(threshold, y, 'k--')

#plt.plot(threshold, echoSBM, 'go', label = r'SBM, 50x20, $p_{cl} = 0.25; p_{add} = 0.0025$, medium mod')
plt.plot(xfSBM, func_lin(xfSBM, *poptSBM), 'g--')
plt.plot(xfSBM1, func_trip(xfSBM1, *poptSBM1), 'g--')
plt.errorbar(threshold, echoSBM, SBMerr, c='g', fmt='o', label = r'SBM, 50x20, $p_{cl} = 0.25; p_{add} = 0.0025$, medium mod')

#plt.plot(threshold, echoSBML, 'yo', label = r'SBM, 100x10, $p_{cl} = 0.06; p_{add} = 0.007$, low mod')
plt.plot(xfSBML, func_lin(xfSBML, *poptSBML), 'y--')
plt.plot(xfSBML1, func_trip(xfSBML1, *poptSBML1), 'y--')
plt.errorbar(threshold, echoSBML, SBMLerr, c='y', fmt='o', label = r'SBM, 100x10, $p_{cl} = 0.06; p_{add} = 0.007$, low mod')'''

#plt.ylim(0, 20)
plt.xlabel("Threshold to change opinion")
plt.ylabel("Fraction of nodes with all neigbors having the same opinion \n(echo chamber size)")
plt.legend(loc='upper right')
plt.title('Size of echo chamber versus threshold to change opinion, 50/50\n N = 1000, REC method \n10 x 10 averaged')
plt.savefig("echo_chamber_vs_threshold_SBM_WS_REC.png")
plt.show()
