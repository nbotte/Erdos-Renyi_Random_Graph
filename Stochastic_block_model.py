import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

block_sizes = [250, 250, 250, 250]
edge_probabilities = [[1., 0.01, 0.01, 0.01], [0.01, 1., 0.01, 0.01], [0.01, 0.01, 1., 0.01], [0.01, 0.01, 0.01, 1.]]
G = nx.generators.community.stochastic_block_model(block_sizes, edge_probabilities, seed=None)
'''degrees = [G.degree(n) for n in G.nodes()]
plt.hist(degrees)
plt.show()'''
print(nx.clustering(G))
print(nx.average_clustering(G))
