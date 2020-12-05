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

class Watts_Strogatz_Network : public Graph{
public:
    // maybe I want the variables private, not sure yet
    int _meanDegree; // the mean degree of the network, should be an even integer (N >> K >> ln(N))
    double _rewireProb; // probability to rewire the edges of the regular lattice

    // constructor
    Watts_Strogatz_Network(int numberOfNodes, int meanDegree, double rewireProb);

    // declare member functions of class Watts-Strogatz
    void makeGraph(); // function that generates a Watts-Strogatz graph
    void makeRegularLattice(); // function that makes a regular lattice
    void rewire(); // function that rewires the edges of the regular lattice    
};

#endif