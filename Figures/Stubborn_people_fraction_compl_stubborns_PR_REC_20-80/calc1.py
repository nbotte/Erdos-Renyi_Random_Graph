import numpy as np

echo = np.loadtxt("Echo_chamber_WS_REC_10-006_fracRes=08_stubb=1_20-80.txt")

begin = echo[:,0]
end = echo[:,1]

print(np.sum(end)/np.sum(begin))
