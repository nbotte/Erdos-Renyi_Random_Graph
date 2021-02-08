import numpy as np
import matplotlib.pyplot as plt

# generate random numbers that sum to fixed value: https://stackoverflow.com/questions/3589214/generate-random-numbers-summing-to-a-predefined-value
# generate power law: https://stackoverflow.com/questions/17882907/python-scipy-stats-powerlaw-negative-exponent/46065079#46065079

# distribution f(x) = C * x^(-a)

# function that generates a power law, k_min is min value for which the power law is generated; k_max the maximum value, gamma is the exponent and y is a random number from the uniform distribution between [0, 1]
def power_law(k_min, k_max, y, gamma):
    return ((k_max**(-gamma+1) - k_min**(-gamma+1))*y  + k_min**(-gamma+1.0))**(1.0/(-gamma + 1.0))

# to generate a distribution you need to create an array
nodes = 1000
scale_free_distribution = np.zeros(nodes)
k_min = 5.0
k_max = 200*k_min
gamma = 1.5

f = open("Community_sizes_powerlaw.txt", "w+")

for n in range(nodes):
    scale_free_distribution[n] = power_law(k_min, k_max,np.random.uniform(0,1), gamma)


hist = np.histogram(scale_free_distribution, bins=1000)
distr = np.sort(scale_free_distribution)
for i in range(len(distr)):
    #print(distr[i])
    distr[i] = int(distr[i])
    #print(distr[i])
'''plt.plot(hist[0])
plt.show()'''

sum = 0
hit = False
sizes = []

while hit is False:
    elem = np.random.choice(distr)
    sum = sum + int(elem)
    if sum == nodes:
        hit = True
        sizes.append(elem)
    elif sum < nodes:
        hit = False
        sizes.append(elem)
    else:
        sum = sum - int(elem)

print(sum)
for elem in sizes:
    f.write(repr(elem) + "\n")
f.close()
