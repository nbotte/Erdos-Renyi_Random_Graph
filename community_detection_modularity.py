import numpy as np
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import networkx.algorithms.community as nx_comm

data = np.loadtxt('facebook_data.txt')

G = nx.Graph()
for i in range(4039):
    G.add_node(i)

edge = np.zeros(2)
for i in range(len(data[:,0])):
    edge[0] = data[:,0][i]
    edge[1] = data[:,1][i]
    G.add_edge(edge[0], edge[1])

c = list(greedy_modularity_communities(G))
mod = nx_comm.modularity(G, c)
print(mod)
