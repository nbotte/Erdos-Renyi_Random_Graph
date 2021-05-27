import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from networkx.generators import random_graphs as gen
from networkx.generators import community as comm


#G = gen.watts_strogatz_graph(50, 4, 0.05)
#G = gen.erdos_renyi_graph(50, 0.08)
n = [10, 10, 10, 10, 10]
p = [[0.3, 0.025, 0.025, 0.025, 0.025], [0.025, 0.3, 0.025, 0.025, 0.025], [0.025, 0.025, 0.3, 0.025, 0.025], [0.025, 0.025, 0.025, 0.3, 0.025], [0.025, 0.025, 0.025, 0.025, 0.3]]
G = comm.stochastic_block_model(n, p)
pos = nx.spring_layout(G)
nx.draw(G, pos, node_size=50)

nx.draw_networkx_nodes(G, pos, nodelist=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], node_color='b', node_size=50)
nx.draw_networkx_nodes(G, pos, nodelist=[10, 11, 12, 13, 14, 15, 16, 17, 18, 19], node_color='g', node_size=50)
nx.draw_networkx_nodes(G, pos, nodelist=[20, 21, 22, 23, 24, 25, 26, 27, 28, 29], node_color='r', node_size=50)
nx.draw_networkx_nodes(G, pos, nodelist=[30, 31, 32, 33, 34, 35, 36, 37, 38, 39], node_color='m', node_size=50)
nx.draw_networkx_nodes(G, pos, nodelist=[40, 41, 42, 43, 44, 45, 46, 47, 48, 49], node_color='y', node_size=50)
plt.draw()
#plt.show(block=False)
#plt.tight_layout()
plt.savefig("SBM.png", dpi=500)
