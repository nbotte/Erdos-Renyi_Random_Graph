import numpy as np

xvalues = np.loadtxt("avDeg_SBM_low.txt")

mean = np.mean(xvalues)
std = np.std(xvalues)
print(mean, std)


'''xAt0op0 = xvalues[:,0]
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

stdAt0op0 = np.std(xAt0op0)
stdAt0op1 = np.std(xAt0op1)
stdAt500op0 = np.std(xAt500op0)
stdAt500op1 = np.std(xAt500op1)

print(varAt0op0, varAt0op1, varAt500op0, varAt500op1)
print(meanAt0op0, meanAt0op1, meanAt500op0, meanAt500op1)
print(stdAt0op0, stdAt0op1, stdAt500op0, stdAt500op1)

X = xAt500op0
X00Float = xAt0op0.astype(float) # convert elements to float (needed for np.reciprocal)

Y = np.reciprocal(X00Float, where=X00Float!=0., dtype=float)

X1 = xAt500op1
X01Float = xAt0op1.astype(float)
Y1 = np.reciprocal(X01Float, where=X01Float!=0., dtype=float)

for i in range(len(Y)):
    if Y[i] > 1000.:
        Y[i] = 0.
    if Y[i] < 0.00000001:
        Y[i] = 0.
for i in range(len(Y1)):
    if Y1[i] > 1000.:
        Y1[i] = 0.
    if Y1[i] < 0.00000001:
        Y1[i] = 0.

E_x = np.mean(X)
E_y = np.mean(Y)

E_x1 = np.mean(X1)
E_y1 = np.mean(Y1)

VarX = np.var(X)
VarY = np.var(Y)

StdX = np.std(X)
StdY = np.std(Y)

VarX1 = np.var(X1)
VarY1 = np.var(Y1)

StdX1 = np.std(X1)
StdY1 = np.std(Y1)

print(E_x, E_y, StdX, StdY)'''

'''for i in range(len(Y)):
    print(X01Float[i])
    print(Y1[i])'''

'''X2 = np.square(X)
Y2 = np.square(Y)

E_X2 = np.mean(X2)
E_Y2 = np.mean(Y2)

X12 = np.square(X1)
Y12 = np.square(Y1)

E_X12 = np.mean(X12)
E_Y12 = np.mean(Y12)

covX2Y2 = np.cov(X2, Y2)[0][1]
covXY = np.cov(X, Y)[0][1]

cov1X2Y2 = np.cov(X12, Y12)[0][1]
cov1XY = np.cov(X1, Y1)[0][1]

VarXY = E_X2 * E_Y2 - E_y**2 * E_x**2

Var1XY = E_X12 * E_Y12 - E_y1**2 * E_x1**2
meanVar = (VarXY + Var1XY)/2
meanSTD = (np.sqrt(VarXY) + np.sqrt(Var1XY))/2
print(np.sqrt(meanVar), meanSTD)
print(E_y*E_x, np.sqrt(VarXY), E_y1*E_x1, np.sqrt(Var1XY))'''
