import numpy as np

xvalues = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_6-001_res=0_xvalues.txt")

xAt0op0 = xvalues[:,0]
xAt0op1 = xvalues[:,1]
xAt500op0 = xvalues[:,2]
xAt500op1 = xvalues[:,3]

varAt0op0 = np.var(xAt0op0)
varAt0op1 = np.var(xAt0op1)
varAt500op0 = np.var(xAt500op0)
varAt500op1 = np.var(xAt500op1)

meanAt0op0 = np.mean(xAt0op0)
meanAt0op1 = np.mean(xAt0op1)
meanAt500op0 = np.mean(xAt500op0)
meanAt500op1 = np.mean(xAt500op1)

print(varAt0op0, varAt0op1, varAt500op0, varAt500op1)

X = xAt500op0
Y = np.reciprocal(xAt0op0)

X1 = xAt500op1
Y1 = np.reciprocal(xAt0op1)

E_x = np.mean(X)
E_y = np.mean(Y)

E_x1 = np.mean(X1)
E_y1 = np.mean(Y1)

VarX = np.var(X)
VarY = np.var(Y)

VarX1 = np.var(X1)
VarY1 = np.var(Y1)

print(E_x, np.mean(xAt0op0), VarX, np.sqrt(VarX), np.var(xAt0op0), np.sqrt(np.var(xAt0op0)))

'''for i in range(len(Y)):
    print(xAt500op0[i])'''

X2 = np.square(X)
Y2 = np.square(Y)

X12 = np.square(X1)
Y12 = np.square(Y1)

covX2Y2 = np.cov(X2, Y2)[0][1]
covXY = np.cov(X, Y)[0][1]

cov1X2Y2 = np.cov(X12, Y12)[0][1]
cov1XY = np.cov(X1, Y1)[0][1]

VarXY = covX2Y2 + (VarX + E_x**2)*(VarY + E_y**2) - (covXY + E_x * E_y)**2

Var1XY = cov1X2Y2 + (VarX1 + E_x1**2)*(VarY1 + E_y1**2) - (cov1XY + E_x1 * E_y1)**2
print(VarXY, np.sqrt(VarXY), Var1XY, np.sqrt(Var1XY))
