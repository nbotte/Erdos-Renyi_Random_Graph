import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.community import greedy_modularity_communities
import networkx.algorithms.community as nx_comm

block_sizes = []
for i in range(int(53)):
    block_sizes.append(151)

edge_probabilities = []
edges = []
for i in range(int(53)):
    for j in range(int(53)):
        if i == j:
            edges.append(0.023)
        else:
            edges.append(0.0001)
    edge_probabilities.append(edges)
    edges = []

G = nx.generators.community.stochastic_block_model(block_sizes, edge_probabilities, seed=None)
'''degrees = [G.degree(n) for n in G.nodes()]
plt.hist(degrees)
plt.show()'''
c = list(greedy_modularity_communities(G))
mod = nx_comm.modularity(G, c)
print(mod)
print(nx.average_clustering(G))
