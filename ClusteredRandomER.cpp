// Nina Botte

#include <cmath>
#include "ClusteredRandomER.h"
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
#include <iterator>
#include <functional>
#include <algorithm>
using namespace std;

// seems to be good!

// for now this class seems to work, but can it be implemented nicer/more efficient? what about use of inheritance?
// TO DO: implement opinion dynamics + make random graphs with number of nodes drawn from log-log distribution

Clustered_Random_Network::Clustered_Random_Network(double rewireProbability){
    _rewireProbability = rewireProbability;
    /* always make sure that you have enough space in your node- and edgelist, otherwise things will need to be copied, but node class has no proper copy constructor, 
    thus will give problems for the neighbours --> can this be done a bit less arbitrary? */
 //   _nodelist.reserve(1500);
   // _edgelist.reserve(10000);
    makeGraph();
}

void Clustered_Random_Network::makeGraph(){
    makeErdosRenyi(250, 0.1, _nodelist.size());
    makeErdosRenyi(250, 0.1, _nodelist.size());
    makeErdosRenyi(250, 0.1, _nodelist.size());
    makeErdosRenyi(250, 0.1, _nodelist.size());

    // rewire edges
    rewireEdges();
}

// function that makes a clustered graph
// maybe add resistance, opinion, active, etc. as arguments so that you can make clusters with different properties
void Clustered_Random_Network::makeErdosRenyi(int numberOfNodes, double edgeProb, int indexStart){
    _numberOfNodes = numberOfNodes;
    _edgeProbability = edgeProb;
    _indexStart = indexStart;
    this->Erdos_Renyi_Network::makeGraph();

  /*  random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    bernoulli_distribution disBern(1.); // make sure that all nodes are active for now

    // add nodes to the graph with some distribution of the 2 possible opinions
    double fractionResistance = 0.; // set the fraction of stubborn/resistant nodes
    double resistance; // variable that determines the resistance of a node
    int opinion; // variable that determines the opinion of a node
    bool active; // variable that determines if node is active
    int index; // variable that determines the index (=name) of a node --> constructed in such a way that the nodes already have the correct name for the complete clustered graph
    for (int i = 0; i < numberOfNodes; i++){
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
        active = disBern(gen);
        index = indexStart + i;
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
    }    
    // add edge between any pair of nodes with a certain probability
    for (int i = 0; i < numberOfNodes; i++){
        for (int j = i+1; j < numberOfNodes; j++){
            double r = dis(gen); // draw a random number that will determine whether there is an edge or not
            if (r < edgeProb){
                int indexIn = indexStart + i;
                int indexOut = indexStart + j;
                Node* N = &_nodelist[indexIn];
                Node* M = &_nodelist[indexOut];
                Edge e = Edge(N, M);
                addEdge(e);
            }
        }
    }*/
}

// function that removes all the edges from the clustered graph
/*void Clustered_Random_Network::removeAllEdges(){
    // remove all the corresponding neighbours
    for (int i = 0; i < _edgelist.size(); i++){
        int indexIn = _edgelist[i].inNode()->index();
        int indexOut = _edgelist[i].outNode()->index();
        _nodelist[indexIn].removeAllNeigh(); // remove neighbours of inNode of edge
        _nodelist[indexOut].removeAllNeigh(); // remove neighbours of outNode of edge
    }
    // Attention: normally clear() would not affect the capacity of the vector, but not always garanteed, so in case of trouble, consider this! (reserve space back after clear)
    _edgelist.clear(); // remove the edges
}

// function that rewires the edges of the different ER graphs with a certain probability
void Clustered_Random_Network::rewireEdges(){
    vector<Edge> helpEdgelist;
    helpEdgelist.reserve(_edgelist.size());

    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    size_t nelem = 1;
    vector<Node> out;
    for (int i = 0; i < _edgelist.size(); i++){
        double r = dis(gen); // random number that will decide if edge is rewired or not
        if (r < _rewireProbability){
            // make sure to compile with c++17 (than sample will not give a problem)
            sample(_nodelist.begin(), _nodelist.end(), back_inserter(out), nelem, mt19937{random_device{}()});
            if (&out.back() != _edgelist[i].inNode()){
                helpEdgelist.push_back(Edge(_edgelist[i].inNode(), &out.back())); // add new edge to the help vector
            }
        }
        else{
            helpEdgelist.push_back(_edgelist[i]);
        }
    }
    removeAllEdges();
    
    for (int i = 0; i < helpEdgelist.size(); i++){
        addEdge(helpEdgelist[i]); // add the new edges in help vector to the edgelist of the clustered graph
    }
}

// function to check if an element is in vector 
bool Clustered_Random_Network::contains(const vector<Edge> vec, Edge e){
    return find(vec.begin(), vec.end(), e) != vec.end();
}

// print function 
void Clustered_Random_Network::print(){
    cout << "Nodes: ";
    for (int i = 0; i < _nodelist.size(); i++){
        cout << _nodelist[i]; 
    }
    cout << endl;
        
    cout << "Edges: ";
    for (int i = 0; i < _edgelist.size(); i++){
        cout << _edgelist[i]; 
    }
    cout << endl;
}*/