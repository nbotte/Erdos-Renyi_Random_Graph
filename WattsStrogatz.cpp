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

// implement default constructor
Watts_Strogatz_Network::Watts_Strogatz_Network(){};

// implement constructor
Watts_Strogatz_Network::Watts_Strogatz_Network(int numberOfNodes, int meanDegree, double rewireProb, double initOp0Frac, int indexStart){
    _numberOfNodes = numberOfNodes;
    _meanDegree = meanDegree;
    _rewireProb = rewireProb;
    _initOp0Frac = initOp0Frac;
    _indexStart = indexStart;

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
    double resistance; // variable that determines the resistance of a node
    int opinion; // variable that determines the opinion of a node
    bool active; // variable that determines if node is active

    vector<int> v;
    // fill a vector with the numbers 0 to numberOfNodes
    for (int i = 0; i < _numberOfNodes; i++){
        v.push_back(i);
    }
    
    // add nodes with a certain random distribution of opinion 0 and 1 to the nodelist
    int N = 0;
    while (v.size()){
        active = 0.; // default: no nodes are active
        resistance = 0.; // default: no stubborn nodes
        
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

    // add edges of the regular lattice
    for (int i = _indexStart; i < _numberOfNodes + _indexStart; i++){
        for (int j = i; j < _numberOfNodes + _indexStart; j++){
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

    for (int i = _indexStart; i < _indexStart + _numberOfNodes; i++){
        for (int index : _nodelist[i].neigh()){
            if (i < index && index <= (i + (_meanDegree/2))){
                double r = dis(gen); // draw a random number that will determine whether the edge is rewired or not
                if (r < _rewireProb){
                    // make sure to compile with c++17 (than sample will not give a problem)
                    // choose a random node that is not equal to i and that is not a neighbor of i
                    sample(_nodelist.begin() + _indexStart, _nodelist.begin() + _indexStart + _numberOfNodes, back_inserter(out), nelem, mt19937{random_device{}()});
                    while (out.back().index() == i || _nodelist[i].containsNeigh(out.back().index())){
                        sample(_nodelist.begin() + _indexStart, _nodelist.begin() + _indexStart + _numberOfNodes, back_inserter(out), nelem, mt19937{random_device{}()});
                    }
                    _nodelist[i].addHelpNeigh(out.back().index());
                }
                else{
                    _nodelist[i].addHelpNeigh(index);
                }
            }
            else if (index > (i + (_meanDegree/2))){
                _nodelist[i].addHelpNeigh(index);
            }
        }
    }
    for (int i = _indexStart; i < _numberOfNodes + _indexStart; i++){
        _nodelist[i].removeAllNeigh();
    }
    ofstream pairFile("Pairs.txt");
    for (int i = _indexStart; i < _numberOfNodes + _indexStart; i++){
        for (int index : _nodelist[i].helpNeigh()){
            _nodelist[i].addNeigh(index);
            pairFile << i << '\t' << index << '\n';
            _nodelist[index].addNeigh(_nodelist[i].index());
        }
        _nodelist[i].removeAllHelpNeigh();
    }
    pairFile.close();
}

// function that calculates the number of triangles in the graph
int Watts_Strogatz_Network::numberOfTriangles(){
    int triangles = 0;
    // loop over all nodes
    for (int i = 0; i < _nodelist.size(); i++){ 
        // for each node: check if there is an edge between any pair of neighbours --> counts each tiangle three times
        for (int index : _nodelist[i].neigh()){
            for (int index2 : _nodelist[i].neigh()){
                // counts each triangle twice
                if (index2 != index){
                    if (checkEdge(_nodelist[index], _nodelist[index2]) || checkEdge(_nodelist[index2], _nodelist[index])){
                        triangles++; 
                    }
                }
            }
        }
    }
    // each triangle is counted three times + each edge is counted twice 
    return triangles/6;
}

// function that calculates the number of connected triples in the graph
int Watts_Strogatz_Network::numberOfTriples(){
    int triples = 0;
    // loop over all nodes
    for (int i  = 0; i < _nodelist.size(); i++){
        // for each node:  count the number of triples with that node in the center
        // https://math.stackexchange.com/questions/60578/what-is-the-term-for-a-factorial-type-operation-but-with-summation-instead-of-p
        int k = _nodelist[i].neigh().size() - 1;
        int triple = k*(k+1)/2;
        triples += triple;
    }
    // each triple is counted twice (the node in the middle of the triplet doesn't count the triplet)
    return triples;
}

// function that returns the overall clustering coefficient (transitivity) of the network
double Watts_Strogatz_Network::overallClustering(){
    double Clus = double(3*numberOfTriangles()) / double(numberOfTriples()); 
    return Clus;
}