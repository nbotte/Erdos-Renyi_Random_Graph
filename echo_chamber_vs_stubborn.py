import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

stubbornSBMH0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_10x100.txt")
stubbornSBMH01 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=01_10x100.txt")
stubbornSBMH02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=02_10x100.txt")
stubbornSBMH03 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=03_10x100.txt")
stubbornSBMH04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=04_10x100.txt")
stubbornSBMH05 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=05_10x100.txt")
stubbornSBMH06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=06_10x100.txt")

stubbornSBM0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=0_50x20.txt")
stubbornSBM01 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=01_50x20.txt")
stubbornSBM02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=02_50x20.txt")
stubbornSBM03 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=03_50x20.txt")
stubbornSBM04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=04_50x20.txt")
stubbornSBM05 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=05_50x20.txt")
stubbornSBM06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_03-0005_res=06_50x20.txt")

stubbornSBML0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=0_100x10.txt")
stubbornSBML01 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=01_100x10.txt")
stubbornSBML02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=02_100x10.txt")
stubbornSBML03 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=03_100x10.txt")
stubbornSBML04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=04_100x10.txt")
stubbornSBML05 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=05_100x10.txt")
stubbornSBML06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-001_res=06_100x10.txt")

stubbornER0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=0.txt")
stubbornER02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=02.txt")
stubbornER04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=04.txt")
stubbornER06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_ER_PR_001_res=06.txt")

stubbornWS0 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_res=0.txt")
stubbornWS01 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_res=01.txt")
stubbornWS02 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_res=02.txt")
stubbornWS03 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_res=03.txt")
stubbornWS04 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_res=04.txt")
stubbornWS05 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_res=05.txt")
stubbornWS06 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_res=06.txt")

def func(x, a, b, c):
    return a*np.exp(b*x) + c

def func_lin(x, a, b):
    return a*x + b

def func_trip(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

echoSBMH = np.zeros(7)

echoSBMH[0] = stubbornSBMH0[:,2][0]
echoSBMH[1] = stubbornSBMH01[:,2][0]
echoSBMH[2] = stubbornSBMH02[:,2][0]
echoSBMH[3] = stubbornSBMH03[:,2][0]
echoSBMH[4] = stubbornSBMH04[:,2][-1]
echoSBMH[5] = stubbornSBMH05[:,2][0]
echoSBMH[6] = stubbornSBMH06[:,2][0]

echoSBM = np.zeros(7)

echoSBM[0] = stubbornSBM0[:,2][0]
echoSBM[1] = stubbornSBM01[:,2][0]
echoSBM[2] = stubbornSBM02[:,2][0]
echoSBM[3] = stubbornSBM03[:,2][0]
echoSBM[4] = stubbornSBM04[:,2][-1]
echoSBM[5] = stubbornSBM05[:,2][-1]
echoSBM[6] = stubbornSBM06[:,2][-1]

echoSBML = np.zeros(7)

echoSBML[0] = stubbornSBML0[:,2][0]
echoSBML[1] = stubbornSBML01[:,2][0]
echoSBML[2] = stubbornSBML02[:,2][0]
echoSBML[3] = stubbornSBML03[:,2][0]
echoSBML[4] = stubbornSBML04[:,2][0]
echoSBML[5] = stubbornSBML05[:,2][0]
echoSBML[6] = stubbornSBML06[:,2][0]

echoER = np.zeros(4)

echoER[0] = stubbornER0[:,2][0]
echoER[1] = stubbornER02[:,2][0]
echoER[2] = stubbornER04[:,2][0]
echoER[3] = stubbornER06[:,2][0]

echoWS = np.zeros(7)

echoWS[0] = stubbornWS0[:,2][-1]
echoWS[1] = stubbornWS01[:,2][-1]
echoWS[2] = stubbornWS02[:,2][-1]
echoWS[3] = stubbornWS03[:,2][-1]
echoWS[4] = stubbornWS04[:,2][-1]
echoWS[5] = stubbornWS05[:,2][-1]
echoWS[6] = stubbornWS06[:,2][-1]

stubborn = np.zeros(7)
for i in range(7):
    stubborn[i] = i/10

poptSBMH, pcovSBMH = curve_fit(func, stubborn, echoSBMH)
print(poptSBMH)
poptSBM, pcovSBM = curve_fit(func_trip, stubborn, echoSBM)
print(poptSBM)
poptSBML, pcovSBML = curve_fit(func, stubborn, echoSBML)
print(poptSBML)
poptER, pcovER = curve_fit(func_lin, stubborn[::2], echoER)
print(poptER)
poptWS, pcovWS = curve_fit(func, stubborn, echoWS)
print(poptWS)
xf = np.linspace(0,0.7,50)

plt.plot(stubborn, echoSBMH, 'o', label = r'SBM, 10x100, $p_{cl} = 0.1; p_{add} = 0.001$, high mod')
plt.plot(xf, func(xf, *poptSBMH))
plt.plot(stubborn, echoSBM, 'o', label = r'SBM, 50x20, $p_{cl} = 0.3; p_{add} = 0.005$, medium mod')
plt.plot(xf, func_trip(xf, *poptSBM))
plt.plot(stubborn, echoSBML, 'o', label = r'SBM, 100x10, $p_{cl} = 0.1; p_{add} = 0.01$, low mod')
plt.plot(xf, func(xf, *poptSBML))
plt.plot(stubborn[::2], echoER, 'o', label = r'ER, $p = 0.01$')
plt.plot(xf, func_lin(xf, *poptER))
plt.plot(stubborn, echoWS, 'o', label = r'WS, $K = 10; \beta = 0.01$')
plt.plot(xf, func(xf, *poptWS))
plt.xlabel("Fraction of stubborn people")
plt.ylabel("Fraction of nodes with all neigbors having the same opinion 0\n(echo chamber size)")
plt.legend(loc='best')
plt.title('Size of echo chamber versus fraction of stubborns, 50/50\n N = 1000, PR method \n10 x 10 averaged')
plt.savefig("echo_chamber_vs_stubborn_SBM_ER_WS_PR.png")
plt.show()
