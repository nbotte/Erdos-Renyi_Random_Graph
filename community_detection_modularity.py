import numpy as np
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import networkx.algorithms.community as nx_comm

data = np.loadtxt('PGP.txt')

G = nx.Graph()
for i in range(10680):
    G.add_node(i)

edge = np.zeros(2)
for i in range(len(data[:,0])):
    edge[0] = data[:,0][i]
    edge[1] = data[:,1][i]
    G.add_edge(edge[0], edge[1])

c = list(greedy_modularity_communities(G))
mod = nx_comm.modularity(G, c)
clus = nx.average_clustering(G)
print(mod)
print(clus)

'''f = open("Communities_real_network.txt", "w+")
for comm in c:
    f.write(repr(comm) + "\n")
f.close()'''
