import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from networkx.generators import random_graphs as gen
from networkx.generators import community as comm
from pyvis.network import Network

net = Network(notebook=True)

#G = gen.watts_strogatz_graph(50, 4, 0.05)
#G = gen.erdos_renyi_graph(50, 0.08)
n = [10, 10, 10, 10, 10]
p = [[0.3, 0.025, 0.025, 0.025, 0.025], [0.025, 0.3, 0.025, 0.025, 0.025], [0.025, 0.025, 0.3, 0.025, 0.025], [0.025, 0.025, 0.025, 0.3, 0.025], [0.025, 0.025, 0.025, 0.025, 0.3]]
G = comm.stochastic_block_model(n, p)
edges = G.edges()
sizes = [15, 15, 15, 15, 15, 15, 15, 15, 15, 15]
net.add_nodes([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], size=sizes, color=['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue'])
net.add_nodes([10, 11, 12, 13, 14, 15, 16, 17, 18, 19], size=sizes, color=['green', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'green'])
net.add_nodes([20, 21, 22, 23, 24, 25, 26, 27, 28, 29], size=sizes, color=['red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red'])
net.add_nodes([30, 31, 32, 33, 34, 35, 36, 37, 38, 39], size=sizes, color=['magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta'])
net.add_nodes([40, 41, 42, 43, 44, 45, 46, 47, 48, 49], size=sizes, color=['cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan'])

net.add_edges(edges)

#net.from_nx(G)

net.show("SBM.html")

#plt.draw()
#plt.show(block=False)
#plt.tight_layout()
#plt.savefig("SBM.png", dpi=500)
