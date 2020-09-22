import numpy as np
from random import random

class Node:
    def __init__(self, index):
        self.index = index # index of node
        self.edges = set() # set of edges attached to node, initially empty

    def addEdge(self, edge):
        self.edges.add(edge) # add edge to edge set

    def removeEdge(self, edge):
        self.edges.discard(edge) # remove edge of edge set

    # function that returns the neighbours of a node
    def neighbours(self):
        return set((edge.inNode if edge.inNode.index != self.index else edge.outNode) for edge in self.edges)

    # print function
    def __str__(self):
        return f"{self.index}"

class Edge:
    def __init__(self, inNode, outNode):
        self.inNode = inNode # start node of edge
        self.outNode = outNode # end node of edge

    # print function
    def __str__(self):
        return f"{self.inNode}-{self.outNode}"

class RandomGraph:
    def __init__(self, numberOfNodes, edgeProbability):
        self.numberOfNodes = numberOfNodes # total number of nodes in graph
        self.edgeProbability = edgeProbability # probability of having an edge between any pair of nodes
        self.nodelist = [] # list of nodes in graph
        self.edgelist = [] # list of edges in graph
        for i in range(numberOfNodes):
            self.addNode(i) # add all the nodes to the graph
        # brute force for now, add edges between nodes with some probability
        for i in range(numberOfNodes):
            for j in range(i, numberOfNodes):
                r = random()
                if r < edgeProbability:
                    node1 = Node(i)
                    node2 = Node(j)
                    self.addEdge(node1, node2)

    # function to add a node to the graph
    def addNode(self, index):
        node = Node(index)
        self.nodelist.append(node)

    # function to add an edge to the graph
    def addEdge(self, inNode, outNode):
        # check if edge is not already there
        if inNode not in outNode.neighbours():
            edge = Edge(inNode, outNode)
            self.edgelist.append(edge)
            inNode.addEdge(edge) # add this new edge to the set of edges attached to inNode
            outNode.addEdge(edge) # add this new edge to the set of edges attached to outNode

    # function to remove an edge from the graph
    def removeEdge(self, edge):
        self.edgelist.remove(edge)
        edge.inNode.removeEdge(edge) # remove edge from the set of edges attached to inNode
        edge.outNode.removeEdge(edge) # remove edge from the set of edges attached to outNode

    # function to remove a node from the graph
    def removeNode(self, node):
        edges = set(node.edges)
        # remove all edges attached the node that needs to be removed
        for edge in edges:
            self.removeEdge(edge)
        self.nodelist.remove(node) # remove node from nodelist

    # print function
    def __str__(self):
        return f"Nodes: {[str(n) for n in self.nodelist]} \n\t Edges: {[str(e) for e in self.edgelist]}"

# Testing
if __name__ == '__main__':
  graph = RandomGraph(100, 0.1)
  node = graph.nodelist[3]
  print(graph)
  print(node.edges)

 # PROBLEM: only prints empty neighbour sets 
