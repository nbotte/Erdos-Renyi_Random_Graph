// Nina Botte

#define _USE_MATH_DEFINES
#include "ErdosRenyi.h"
#include "Graph.h"
#include <cmath>
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <string>
#include <vector>
using namespace std;

#ifndef CLUSTEREDRANDOMER_H
#define CLUSTEREDRANDOMER_H

class Clustered_Random_Network : public Erdos_Renyi_Network{
    double _rewireAddProbability; // probability of rewireing or adding edges between clusters (whether you rewire or add depends on 'type' parameter in constructor)
    string _type; // tells whether you rewire or you add edges to make the clustered graph 

public:
    // define a constructor
    Clustered_Random_Network(double rewireAddProbability, string type);

    // declare member functions of class Clustered_Random_Network 
    vector<int> makeErdosRenyi(int numberOfNodes, double edgeProb, int indexStart); // function that makes an ER-graph + returns a vector with the indices of the nodes in that cluster
    void makeGraph(); // function that makes a clustered graph
    void rewireEdges(vector<vector<int>>); // function that rewires the edges of a graph, takes a vector of vectors of indices of the nodes of the different clusters as argument
    void addEdges(vector<vector<int>>);
};

#endif