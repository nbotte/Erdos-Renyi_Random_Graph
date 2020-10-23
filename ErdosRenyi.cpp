// Nina Botte

#include <cmath>
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

// constructor, construct graph by making the nodes and the edges with a given probability
Erdos_Renyi_Network::Erdos_Renyi_Network(int numberOfNodes, double edgeProbability, double bernouilliProbability, int indexStart){
    _numberOfNodes = numberOfNodes;
    _edgeProbability = edgeProbability;
    _bernouilliProbability = bernouilliProbability;

    // reserve enough memory space for the vectors
    _nodelist.reserve(_numberOfNodes); 
    _edgelist.reserve(pow(_numberOfNodes, 2)); 

    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    bernoulli_distribution disBern(_bernouilliProbability);

    // add nodes to the graph with some distribution of the 2 possible opinions
    double fractionResistance = 0.; // set the fraction of stubborn/resistant nodes
    double resistance; // variable that determines the resistance of a node
    int opinion; // variable that determines the opinion of a node
    bool active; // variable that determines if node is active
    int index; // variable that gives the index (name) of a node
    for (int i = 0; i < _numberOfNodes; i++){
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
        index = indexStart + i; // indexStart is normally zero, but not if we need to construct clustered graphs out of different Erdos-Renyi graphs
        Node n = Node(index, opinion, resistance, active);
        addNode(n);
    }    
    // add edge between any pair of nodes with a certain probability
    for (int i = 0; i < _nodelist.size(); i++){
        for (int j = i+1; j < _nodelist.size(); j++){
            double r = dis(gen); // draw a random number that will determine whether there is an edge or not
            if (r < _edgeProbability){
                int indexIn = i;
                int indexOut = j;
                addEdge(indexIn, indexOut);
            }
        }
    }
}

// getter, provides access to data member with corresponding name
vector<Node> Erdos_Renyi_Network::nodelist() {return _nodelist;}
vector<Edge> Erdos_Renyi_Network::edgelist() {return _edgelist;}

// function to add a node to the graph
void Erdos_Renyi_Network::addNode(Node n){
    _nodelist.push_back(n);
}

// function to add an edge to the graph
void Erdos_Renyi_Network::addEdge(int indexIn, int indexOut){
    Node* N = &_nodelist[indexIn];
    Node* M = &_nodelist[indexOut];
    Edge e = Edge(N, M);
    // check if edge is already there
    if (contains(_edgelist, e) == false){
        _edgelist.push_back(e);
        _nodelist[indexIn].addNeigh(&_nodelist[indexOut]); // add outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].addNeigh(&_nodelist[indexIn]); // add inNode of edge to neighbours of outNode of edge

    }
}

// function to change the opinions of the nodes in graph based on majority model
void Erdos_Renyi_Network::changeOpinions(){
    // give all the nodes a new opinion
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].changeOpinion();
    }
    // set the opinion of the nodes equal to their new opinion
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].setNewOpinion();
    }
}

// function to deactivate all the nodes in the network, but first the current active nodes need to set their wasActive variable to true
// this function might be unneccessary!
void Erdos_Renyi_Network::deactivateNodes(){
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].setWasActive();
        _nodelist[i].deactivate();
    }
}

// function that sets a fraction of the nodes as active (according to a bernouilli distribution)
void Erdos_Renyi_Network::setNodesActive(){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    bernoulli_distribution disBern(_bernouilliProbability);
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].setActive(disBern(gen));
    }
}

// function to remove node and edge are not good yet, but are they needed?
// function to remove an edge from the graph
// NOT TESTED YET
/*void removeEdge(Edge* e){
    _edgelist.erase(remove(_edgelist.begin(), _edgelist.end(), e), _edgelist.end()); // remove edge from edgelist
    e->inNode()->removeNeigh(e->outNode()->index()); // remove outNode from the set of neighbours attached to inNode
    e->outNode()->removeNeigh(e->inNode()->index()); // remove inNode from the set of neighbours attached to outNode
}

// function to remove a node from the graph
// NOT TESTED YET
void removeNode(int index){
    // remove all neighbours from the neighbour list of that node
    _nodelist[index].neigh().clear();
    // remove node from nodelist
    _nodelist.erase(index);
}*/

// function to check if an element is in vector 
bool Erdos_Renyi_Network::contains(const vector<Edge> vec, Edge e){
    return find(vec.begin(), vec.end(), e) != vec.end();
}

// function that counts the fraction of the 2 opinions in the graph (returns a vector with the 2 fractions)
vector<double> Erdos_Renyi_Network::countOpinionFraction(){
    int opinion0 = 0;
    int opinion1 = 0;
    vector<double> fractions;
    for (int i = 0; i < _nodelist.size(); i++){
        if (_nodelist[i].opinion() == 0){
            opinion0++;
        }
        else if (_nodelist[i].opinion() == 1){
            opinion1++;
        }
    }
    fractions.push_back(double(opinion0)/double(_numberOfNodes));
    fractions.push_back(double(opinion1)/double(_numberOfNodes));
    return fractions;
}

// print function 
void Erdos_Renyi_Network::print(){
    cout << "Nodes: ";
    for (int i = 0; i < _nodelist.size(); ++i){
        cout << _nodelist[i]; 
    }
    cout << endl;
        
    cout << "Edges: ";
    for (int i = 0; i < _edgelist.size(); i++){
        cout << _edgelist[i]; 
    }
    cout << endl;
}

