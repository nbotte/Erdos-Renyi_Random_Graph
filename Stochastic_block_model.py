import numpy as np
import networkx as nx

block_sizes = [25, 25, 25, 25]
edge_probabilities = [[0.1, 0.01, 0.01, 0.01], [0.01, 0.1, 0.01, 0.01], [0.01, 0.01, 0.1, 0.01], [0.01, 0.01, 0.01, 0.1]]
G = nx.generators.community.stochastic_block_model(block_sizes, edge_probabilities, seed=None)

print(G.edges())
