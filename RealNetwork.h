// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

#define _USE_MATH_DEFINES
#include <cmath>
#include "Graph.h"
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <vector>
using namespace std;

#ifndef REALNETWORK_H
#define REALNETWORK_H

class Real_World_Network : virtual public Graph{
public:
    vector<vector<int>> _edges; // vector that contains the edges (= vector with the inNode and the outNode of each edge)

    // constructor
    Real_World_Network(int numberOfNodes, vector<vector<int>> edges);

    // declare member functions of class Real_World_Network
    void makeGraph(); // function that generates the corresponding real-world graph
    double calculateModularity(vector<vector<int>> communities); // function that calculates the modularity (returns the calculated modularity)
};

#endif