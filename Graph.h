// Nina Botte

#define _USE_MATH_DEFINES
#include <cmath>
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <vector>
using namespace std;

#ifndef GRAPH_H
#define GRAPH_H

// Abstract base class: Graph

class Graph{
public: 
    // define a virtual destructor
    virtual ~Graph();

    vector<Node> _nodelist; // vector of nodes in graph
    vector<Edge> _edgelist; // vector of edges in graph
    int _numberOfNodes; // total number of nodes in the graph
    double _rewireProbability; // probability of rewireing edges 

    // define getters, provides access to data member with corresponding name
    vector<Node> nodelist();
    vector<Edge> edgelist();

    // define a virtual function makeGraph
    virtual void makeGraph() = 0;

    // declare member function of graph
    void addNode(Node n); // function to add nodes to the nodelist
    void addEdge(Edge e); // function to add edges to the edgelist
    void removeEdge(Edge e); // function to remove edges from the edgelist
    void removeAllEdges(); // function that removes all the edges from the graph
    void rewireEdges(); // function that rewires the edges of a graph

    // opinion dynamics functions of a graph
    void changeOpinions();
    void deactivateNodes();
    void setNodesActive();
    vector<double> countOpinionFraction();

    // print function for a graph
    void print();

    // not a member function, will be used to check if an edge is already in the edgelist --> should be written outside class definition?
    bool contains(const vector<Edge> vec, Edge e);

};


#endif