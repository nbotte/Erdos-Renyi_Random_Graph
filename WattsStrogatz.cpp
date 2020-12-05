// Nina Botte

#include <cmath>
#include <memory>
#include "Graph.h"
#include "Node.h"
#include "Edge.h"
#include "WattsStrogatz.h"
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

// https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model

// implement constructor
Watts_Strogatz_Network::Watts_Strogatz_Network(int numberOfNodes, int meanDegree, double rewireProb){
    _numberOfNodes = numberOfNodes;
    _meanDegree = meanDegree;
    _rewireProb = rewireProb;

    _nodelist.resize(_numberOfNodes); 

    // make Watts-Strogatz network by calling the function makeGraph()
    makeGraph();
}

// function that makes a Watts-Strogatz network
void Watts_Strogatz_Network::makeGraph(){
    makeRegularLattice();
    rewire();
}

// function that makes a regular lattice
void Watts_Strogatz_Network::makeRegularLattice(){
    // opinion, active, resistance = 0 for now
    double resistance = 0.; // variable that determines the resistance of a node
    int opinion = 0; // variable that determines the opinion of a node
    bool active = false; // variable that determines if node is active
    
    // add nodes to the nodelist
    for (int i = 0; i < _numberOfNodes; i++){
        int index = i;
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
    }

    // add edges of the regular lattice
    for (int i = 0; i < _numberOfNodes; i++){
        for (int j = i; j < _numberOfNodes; j++){
           if (0 < abs(i-j)%(_numberOfNodes-1-(_meanDegree/2)) && abs(i-j)%(_numberOfNodes-1-(_meanDegree/2)) <= (_meanDegree/2)){
                auto N = make_shared<Node>(_nodelist[i]);
                auto M = make_shared<Node>(_nodelist[j]);
                Edge e = Edge(N, M);
                addEdge(e);
           } 
        }
    }
}

// function that rewires the edges of the regular lattice
void Watts_Strogatz_Network::rewire(){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    size_t nelem = 1;
    vector<Node> out; // will contain the random node to which you rewire

    for (int i = 0; i < _nodelist.size(); i++){
        for (int index : _nodelist[i].neigh()){
            if (i < index && index <= (i + (_meanDegree/2))){
                double r = dis(gen); // draw a random number that will determine whether the edge is rewired or not
                if (r < _rewireProb){
                    // make sure to compile with c++17 (than sample will not give a problem)
                    // choose a random node that is not equal to i and that is not a neighbor of i
                    sample(_nodelist.begin(), _nodelist.end(), back_inserter(out), nelem, mt19937{random_device{}()});
                    while (out.back().index() == i || _nodelist[i].containsNeigh(out.back().index())){
                        sample(_nodelist.begin(), _nodelist.end(), back_inserter(out), nelem, mt19937{random_device{}()});
                    }
                    _nodelist[i].removeNeigh(index);
                    _nodelist[i].addNeigh(out.back().index());
                }
            }
        }
    }
}