// Nina Botte

#define _USE_MATH_DEFINES
#include <cmath>
#include "Graph.h"
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <vector>
using namespace std;

#ifndef WATTSSTROGATZ_H
#define WATTSSTROGATZ_H

class Watts_Strogatz_Network : virtual public Graph{
public:
    // maybe I want the variables private, not sure yet
    int _meanDegree; // the mean degree of the network, should be an even integer (N >> K >> ln(N))
    double _rewireProb; // probability to rewire the edges of the regular lattice
    double _initOp0Frac; // initial fraction of nodes with opinion 0
    int _indexStart; // this is mainly needed to construct WS networks as constituents for the SBM_WS

    // define default constructor
    Watts_Strogatz_Network();

    // constructor
    Watts_Strogatz_Network(int numberOfNodes, int meanDegree, double rewireProb, double initOp0Frac, int indexStart);

    // declare member functions of class Watts-Strogatz
    void makeGraph(); // function that generates a Watts-Strogatz graph
    void makeRegularLattice(); // function that makes a regular lattice
    void rewire(); // function that rewires the edges of the regular lattice  

    int numberOfTriangles(); // function that calculates the number of triangles in the graph
    int numberOfTriples(); // function that calculates the number of connected triples in the graph

    double overallClustering(); // function that returns the overall clustering coefficient (transitivity) of the network

};

#endif