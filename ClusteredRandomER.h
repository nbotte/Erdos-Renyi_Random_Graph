// Nina Botte

#define _USE_MATH_DEFINES
#include <cmath>
#include "Node.h"
#include "Edge.h"
#include "ErdosRenyi.h"
#include <math.h>
#include <vector>
using namespace std;

#ifndef CLUSTEREDRANDOMER_H
#define CLUSTEREDRANDOMER_H

class Clustered_Random_Network{
    double _rewireProbability; // probability of rewireing the edges of the different ER-graphs that make up the clustered graph
    vector<Node> _nodelist; // vector of nodes in clustered graph
    vector<Edge> _edgelist; // vector of edges in clustered graph

public:
    // define a constructor
    Clustered_Random_Network(double rewireProbability);

    // define getters, provides access to data member with corresponding name
    vector<Node> nodelist() const;
    vector<Edge> edgelist() const;

    // declare member functions of class Clustered_Random_Network 
    // opinion dynamic functions not included yet! Or maybe make a different class for the opinion dynamics stuff?
    void addNode(Node n); // function to add nodes to the nodelist
    void addEdge(Edge e); // function to add edges to the edgelist
    void removeEdge(Edge e); // function to remove edges from the edgelist
    void removeAllEdges(); // function to remove all the edges from the graph
    void rewireEdges();
    void print();

    // not a member function, will be used to check if an edge is already in the edgelist
    bool contains(const vector<Edge> vec, Edge e);
};

#endif