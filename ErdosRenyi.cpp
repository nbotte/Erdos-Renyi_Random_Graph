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

    // add nodes to the graph with some distribution of the 2 possible opinions
    double fractionResistance = 0.; // set the fraction of stubborn/resistant nodes
    double resistance; // variable that determines the resistance of a node
    int opinion; // variable that determines the opinion of a node
    bool active; // variable that determines if node is active

    vector<int> v;
    // fill a vector with the numbers 0 to numberOfNodes
    for (int i = 0; i < _numberOfNodes; i++){
        v.push_back(i);
    }

    /*for (int i = 0; i < _numberOfNodes; i++){
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
    }*/

    // first 500 indices obtained from the permutated vector v will have an opinion zero, others get opinion 1
    // resistance of node is not implemented yet
    int N = 0;
    while (v.size()){
        active = 0.; // default: no nodes are active
        resistance = 0.; // good for now, no stubborn nodes
        int index = getRandomElement(v, _numberOfNodes - 1) + _indexStart;
        // Note: _numberOfNodes should be even, otherwise you will get a bias!
        if (N < _numberOfNodes/2){
            opinion = 0;
        }
        else{
            opinion = 1;
        }
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
        N++;
    }
   /* for (int i = 0; i < _numberOfNodes/2; i++){
        double k = dis(gen); // random number to determine if node is stubborn
        if (k <= fractionResistance){
            resistance = 0.;
        }
        else{
            resistance = 0.;
        }
        active = 0.;
        int index = _indexStart + i;
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
    } */
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

int getRandomElement(vector<int>& v, int length){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_int_distribution<> dis(0, length); // uniform diwtribution between 0 and numberOfNodes (here: 1000)
    int n = v.size();
    int index = dis(gen) % n; // random number between 0 and 999 --> but make sure that it is always in the range of v (size of v changes!)
    int elem = v[index]; // get random element from vector
    swap(v[index], v[n-1]);
    v.pop_back();
    return elem;    
}



