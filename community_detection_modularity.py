import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import networkx.algorithms.community as nx_comm

data = np.loadtxt('lastfm_fin_edgelist.txt')

G = nx.Graph()
for i in range(8003):
    G.add_node(i)

edge = np.zeros(2)
for i in range(len(data[:,0])):
    edge[0] = data[:,0][i]
    edge[1] = data[:,1][i]
    G.add_edge(edge[0], edge[1])

c = list(greedy_modularity_communities(G))
'''mod = nx_comm.modularity(G, c)
clus = nx.average_clustering(G)
print(mod)
print(clus)'''
sizes = []
for comm in c:
    sizes.append(len(comm))

f = open("Communities_sizes_real_network_lastfm.txt", "w+")
for size in sizes:
    f.write(repr(size) + "\n")
f.close()

'''bins = np.arange(0, 1800, 1)

#plt.xlim([min(sizes)-5, max(sizes)+5])

plt.hist(sizes, bins=bins)
#plt.title('Random Gaussian data (fixed bin size)')
plt.xlabel('variable X (bin size = 5)')
plt.ylabel('count')

plt.show()'''

'''f = open("Communities_real_network.txt", "w+")
for comm in c:
    f.write(repr(comm) + "\n")
f.close()'''
