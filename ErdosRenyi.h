// Nina Botte

#define _USE_MATH_DEFINES
#include <cmath>
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <vector>
using namespace std;

#ifndef ERDOSRENYI_H
#define ERDOSRENYI_H

class Erdos_Renyi_Network{
    int _numberOfNodes; // total number of nodes in the graph
    double _edgeProbability; // the probablility of having an edge between any pair of nodes
    vector<Node> _nodelist; // vector of nodes in graph
    vector<Edge> _edgelist; // vector of edges in graph
    double _bernouilliProbability; // probability of having true in bernouilli process (will determine if a node is active or not)

public:
    // define a constructor
    Erdos_Renyi_Network(int, double, double);

    // define getters, provides access to data member with corresponding name
    vector<Node> nodelist();
    vector<Edge> edgelist();

    // declare member functions of class Erdos-Renyi
    void addNode(Node n);
    void addEdge(int indexIn, int indexOut);
    void changeOpinions();
    void deactivateNodes();
    void setNodesActive();
    vector<double> countOpinionFraction();
    void print();

    // not a member function
    bool contains(const vector<Edge> vec, Edge e);

};


#endif