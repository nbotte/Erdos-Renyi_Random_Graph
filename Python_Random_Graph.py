import numpy as np
from random import random

class Node:
    def __init__(self, index):
        self.index = index # index of node
        self.edges = [] # list of edges attached to node, initially empty

    def addEdge(self, edge):
        self.edges.append(edge) # add list to edge list

    def removeEdge(self, edge):
        self.edges.remove(edge) # remove list to edge list

class Edge:
    def __init__(self, inNode, outNode):
        self.inNode = inNode # start node of edge
        self.outNode = outNode # end node of edge

class RandomGraph:
    def __init__(self, numberOfNodes, edgeProbability):
        self.numberOfNodes = numberOfNodes # total number of nodes in graph
        self.edgeProbability = edgeProbability # probability of having an edge between any pair of nodes
        self.nodelist = []
        self.edgelist = []
        for i in range(numberOfNodes):
            self.addNode(i)
        # brute force for now
        for i in range(numberOfNodes):
            for j in range(numberOfNodes):
                if i != j:
                    r = random()
                    if r < edgeProbability:
                        node1 = Node(i)
                        node2 = Node(j)
                        self.addEdge(node1, node2)


    def getNodelist(self):
        return self.nodelist

    def getEdgelist(self):
        return self.edgelist

    def addNode(self, index):
        node = Node(index)
        self.nodelist.append(node)

    def addEdge(self, node1, node2):
        edge = Edge(node1, node2)
        self.edgelist.append(edge)
        node1.addEdge(edge)
        node2.addEdge(edge)

    def removeEdge(self, edge):
        self.edgelist.remove(edge)
        inNode.removeEdge(edge)
        outNode.remove(edge)

    def removeNode(self, node):
        edges = node.edges()
        for edge in edges:
            self.removeEdge(edge)
        self.nodelist.remove(node)

# Testing
if __name__ == '__main__':
  graph = RandomGraph(100, 0.5)
  print(graph.getEdgelist())
