// Nina Botte

#include <cmath>
#include <memory>
#include "Graph.h"
#include "ErdosRenyi.h"
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <random>
#include <functional>
#include <algorithm>
using namespace std;

// implement default constructor
Erdos_Renyi_Network::Erdos_Renyi_Network(){};

// constructor, construct graph by making the nodes and the edges with a given probability
Erdos_Renyi_Network::Erdos_Renyi_Network(int numberOfNodes, double edgeProbability, double bernouilliProbability, int indexStart){
    _numberOfNodes = numberOfNodes;
    _edgeProbability = edgeProbability;
    _bernouilliProbability = bernouilliProbability;
    _indexStart = indexStart;

    // reserve enough memory space for the vectors
    _nodelist.reserve(_numberOfNodes); 
    _edgelist.reserve(pow(_numberOfNodes, 2)); 

    // make random graph by calling the function makeGraph
    makeGraph();    
}

void Erdos_Renyi_Network::makeGraph(){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    bernoulli_distribution disBern(_bernouilliProbability);

    // add nodes to the graph with some distribution of the 2 possible opinions
    double fractionResistance = 0.; // set the fraction of stubborn/resistant nodes
    double resistance; // variable that determines the resistance of a node
    int opinion; // variable that determines the opinion of a node
    bool active; // variable that determines if node is active
    
    for (int i = 0; i < _numberOfNodes; i++){
        if (i % 2){
            opinion = 0;
        }
        else{
            opinion = 1;
        }
        resistance = 0.;
        active = 0.;
        int index = _indexStart + i;
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
    }

  /*  for (int i = 0; i < _numberOfNodes; i++){
        double k = dis(gen); // random number to determine if node is stubborn
        if (k <= fractionResistance){
            resistance = 0.;
        }
        else{
            resistance = 0.;
        }
        double r = dis(gen); // random number to determine the opinion of a node
        if (r < 0.5){
            opinion = 0;
        }
        else{
            opinion = 1;
        }
        active = 0.;
        int index = _indexStart + i;
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
    } */
    // add edge between any pair of nodes with a certain probability
    for (int i = 0; i < _nodelist.size() - _indexStart; i++){
        for (int j = i+1; j < _nodelist.size() - _indexStart; j++){
            double r = dis(gen); // draw a random number that will determine whether there is an edge or not
            if (r < _edgeProbability){
                int indexIn = _indexStart + i;
                int indexOut = _indexStart + j;
                auto N = make_shared<Node>(_nodelist[indexIn]);
                auto M = make_shared<Node>(_nodelist[indexOut]);
                Edge e = Edge(N, M);
                addEdge(e);
            }
        }
    }
}



