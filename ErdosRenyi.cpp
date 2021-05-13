// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

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
Erdos_Renyi_Network::Erdos_Renyi_Network(int numberOfNodes, double edgeProbability, double bernouilliProbability, double initOp0Frac, int indexStart){
    _numberOfNodes = numberOfNodes; // set the total number of nodes in the graph
    _edgeProbability = edgeProbability; // set the edge probability
    _bernouilliProbability = bernouilliProbability; // set the bernouilli probability for being active
    _indexStart = indexStart; // set the index from which the node names start --> mainly used when constructing a SBM where each community is made by a seperate ER graph
    _initOp0Frac = initOp0Frac; // set the initial fraction of nodes with opinion 0

    // reserve enough memory space for the vectors
    _nodelist.resize(_numberOfNodes); 
    _edgelist.reserve(pow(_numberOfNodes, 2)); // not needed?

    // make random graph by calling the function makeGraph
    makeGraph();    
}

void Erdos_Renyi_Network::makeGraph(){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    bernoulli_distribution disBern(_bernouilliProbability);

    double resistance; // variable that determines the resistance of a node
    int opinion; // variable that determines the opinion of a node
    bool active; // variable that determines if node is active

    vector<int> v;
    // fill a vector with the numbers 0 to numberOfNodes
    for (int i = 0; i < _numberOfNodes; i++){
        v.push_back(i);
    }

    // first 500 indices obtained from the permutated vector v will have an opinion zero, others get opinion 1
    int N = 0;
    while (v.size()){
        active = 0.; // default: no nodes are active
        double r = dis(gen); // random number that will determine if node is stubborn or not
        
        resistance = 0.; // set default resistance to zero
               
        int index = getRandomElement(v, _numberOfNodes - 1) + _indexStart;
        // Note: _numberOfNodes should be even, otherwise you will get a bias!
        if (N < int(_numberOfNodes*_initOp0Frac)){
            opinion = 0;
        }
        else{
            opinion = 1;
        }
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
        N++;
    }

    // add edge between any pair of nodes with a certain probability
    for (int i = 0; i < _numberOfNodes; i++){
        for (int j = i+1; j < _numberOfNodes; j++){
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





