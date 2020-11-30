// Nina Botte

#define _USE_MATH_DEFINES
#include <cmath>
#include "Graph.h"
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <vector>
using namespace std;

#ifndef ERDOSRENYI_H
#define ERDOSRENYI_H

class Erdos_Renyi_Network : public Graph{
public:

    double _edgeProbability; // the probablility of having an edge between any pair of nodes
    double _bernouilliProbability; // probability of having true in bernouilli process (will determine if a node is active or not)
    int _indexStart; 

    // define default constructor
    Erdos_Renyi_Network();

    // define a constructor
    Erdos_Renyi_Network(int, double, double, int);

    // declare member functions of class Erdos-Renyi
    void makeGraph(); // function that makes an Erdos-Renyi graph

};

// function than returns a random element from a vector
int getRandomElement(vector<int>& v, int length);


#endif