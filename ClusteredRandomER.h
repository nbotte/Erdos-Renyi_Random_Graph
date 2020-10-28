// Nina Botte

#define _USE_MATH_DEFINES
#include "ErdosRenyi.h"
#include "Graph.h"
#include <cmath>
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <vector>
using namespace std;

#ifndef CLUSTEREDRANDOMER_H
#define CLUSTEREDRANDOMER_H

class Clustered_Random_Network : public Erdos_Renyi_Network{

public:
    // define a constructor
    Clustered_Random_Network(double rewireProbability);

    // declare member functions of class Clustered_Random_Network 
    void makeErdosRenyi(int numberOfNodes, double edgeProb, int indexStart); // function that makes an ER-graph
    void makeGraph(); // function that makes a clustered graph
};

#endif