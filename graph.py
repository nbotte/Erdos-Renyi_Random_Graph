import numpy as np


class Edge:
    def __init__(self, node_out, node_in):
        self.node_out = node_out
        self.node_in = node_in

    def __str__(self):
        return f"{self.node_out}-{self.node_in}"


class Node:
    def __init__(self, name):
        self.name = name
        self.outgoing_edges = set()

    def add_edge(self, edge):
        self.outgoing_edges.add(edge)

    def remove_edge(self, edge):
        self.outgoing_edges.discard(edge)

    @property
    def neighbours(self):
        return set((edge.node_in if edge.node_in.name != self.name else edge.node_out) for edge in self.outgoing_edges)

    def __str__(self):
        return f"({self.name})"


class UndirectedGraph:
    def __init__(self, node_names, adjacency_matrix=None):
        self._edges = []
        self._nodes = [Node(name) for name in node_names]
        if adjacency_matrix is not None:
            length = len(adjacency_matrix)
            for i in range(length):
                for j in range(i, length):
                    if adjacency_matrix[i, j]:
                        new_edge = Edge(self._nodes[i], self._nodes[j])
                        self._edges.append(new_edge)
                        self._nodes[i].add_edge(new_edge)
                        self._nodes[j].add_edge(new_edge)

    def add_node(self, name):
        self._nodes.append(Node(name))

    def add_edge(self, node_out, node_in):
        if node_in not in node_out.neighbours:
            new_edge = Edge(node_out, node_in)
            node_out.add_edge(new_edge)
            node_in.add_edge(new_edge)
            self._edges.append(new_edge)
            print("added edge")
        else:
            print("edge already present")

    def remove_node(self, node):
        for edge in set(node.outgoing_edges):
            self.remove_edge(edge)
        self._nodes.remove(node)

    def remove_edge(self, edge):
        self._edges.remove(edge)
        edge.node_in.remove_edge(edge)
        edge.node_out.remove_edge(edge)

    def __str__(self):
        return f"UndirectedGraph with\n\tNodes: {[str(n) for n in self._nodes]}\n\tEdges: " \
               f"{[str(e) for e in self._edges]}"

    @property
    def adjacency_matrix(self):
        length = len(self._nodes)
        adjacency_matrix = np.zeros((length, length), np.int0)
        for i in range(length):
            for j in range(i, length):
                if self._nodes[j] in self._nodes[i].neighbours:
                    adjacency_matrix[i, j] = 1
                    adjacency_matrix[j, i] = 1
        return adjacency_matrix


if __name__ == "__main__":
    np.random.seed(0)
    am = np.random.random((5, 5))
    am[np.diag_indices(5)] = 0
    am = ((am + am.T) / 2 > 0.5).astype(np.int0)
    print(am)
    print("")
    names = [f"node_{i+1}" for i in range(5)]
    graph = UndirectedGraph(node_names=names, adjacency_matrix=am)
    print(graph)

    n = graph._nodes[0]
    print(f"\ndeleting node: {n}")
    graph.remove_node(n)
    print(graph)

    print("\nAdding node \"bla\"")
    graph.add_node("bla")
    print(graph)

    n1 = graph._nodes[3]
    n2 = graph._nodes[-1]
    print(f"\nAdding edge from node {n1} to node {n2}")
    graph.add_edge(n1, n2)
    print(graph)

    print(f"\nAdding edge from node {n1} to node {n2}")
    graph.add_edge(n1, n2)
    print(graph)

    e = graph._edges[1]
    print(f"\nDeleting edge: {e}")
    graph.remove_edge(e)
    print(graph)

    print(f"\nAdjacency matrix after modifications is: \n{graph.adjacency_matrix}")
