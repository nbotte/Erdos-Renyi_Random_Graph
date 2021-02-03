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


    int _numberOfNodes; // total number of nodes in the graph
    vector<Node> _nodelist; // vector of nodes in graph
    vector<Edge> _edgelist; // vector of edges in graph --> Not needed anymore?

    // define getters, provides access to data member with corresponding name
    vector<Node> nodelist();
    vector<Edge> edgelist();

    // define a virtual function makeGraph --> each subclass of Graph will have a different implementation of this function
    virtual void makeGraph() = 0;

    // declare member function of graph
    void addNode(Node n); // function to add nodes to the nodelist
    void addEdge(Edge e); // function to add edges to the edgelist
    void removeEdge(Edge e); // function to remove edges from the edgelist
    void removeAllEdges(); // function that removes all the edges from the graph
    bool checkEdge(Node u, Node v); // function that checks whether there is an edge between 2 vertices u and v that are both neigbors of the same node (used to calculate the local clustering coefficient), returns True if there is an edge
    double localClustering(Node u); // function that calculates the local clustering coefficient of a node in the graph
    double averageClustering(); // function that calculates the average clustering coefficient of the graph
    int numberOfEdges(); // function that calculates the number of edges in the graph

    // opinion dynamics functions of a graph
    void changeOpinions(); // changes the opinions of the nodes in the graph
    void deactivateNodes(); // makes all the nodes inactive
    void setNodesActive(double bernProb); // sets nodes active according to a bernouillidistribution
    void resetInitOpinion(double initOp0Frac); // resets the initial opinions of the nodes in the graph
    vector<double> countOpinionFraction();

    void changeRandomOpinion();

    // print function for a graph
    void print();

    // not a member function, will be used to check if an edge is already in the edgelist --> should be written outside class definition?
    bool contains(const vector<Edge> vec, Edge e);
};

// function than returns a random element from a vector
int getRandomElement(vector<int>& v, int length);

#endif