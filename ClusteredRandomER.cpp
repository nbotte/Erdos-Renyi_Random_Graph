// Nina Botte

#include <cmath>
#include "ClusteredRandomER.h"
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

// for now this class seems to work, but can it be implemented nicer/more efficient? what about use of inheritance?

// No opinion dynamics implemented yet --> to do

Clustered_Random_Network::Clustered_Random_Network(double rewireProbability){
    _rewireProbability = rewireProbability;
    /* always make sure that you have enough space in your node- and edgelist, otherwise things will need to be copied, but node class has no proper copy constructor, 
    thus will give problems for the neighbours */
    _nodelist.reserve(150);
    _edgelist.reserve(200);

    cout << _nodelist.size();
    makeErdosRenyi(25, 0.1, _nodelist.size());
    cout << _nodelist.size();
    makeErdosRenyi(25, 0.1, _nodelist.size());
    cout << _nodelist.size();
    makeErdosRenyi(25, 0.1, _nodelist.size());
    cout << _nodelist.size();
    makeErdosRenyi(25, 0.1, _nodelist.size());

    // rewire edges
    rewireEdges();
}

// implementation of the getters
vector<Node> Clustered_Random_Network::nodelist() const {return _nodelist;}
vector<Edge> Clustered_Random_Network::edgelist() const {return _edgelist;}

// function that adds a node to the nodelist of the clustered graph NEEDED?
void Clustered_Random_Network::addNode(Node n){
    _nodelist.push_back(n);
}

// function that adds an edge to the edgelist of the clustered graph + adds the corresponding neighbours to the involved nodes
void Clustered_Random_Network::addEdge(Edge e){
    // check if edge is already there
    if (contains(_edgelist, e) == false){
        _edgelist.push_back(e);
        int indexIn = e.inNode()->index();
        int indexOut = e.outNode()->index();
        //cout << indexIn << ' ' << indexOut << endl;
        _nodelist[indexIn].addNeigh(&_nodelist[indexOut]); // add outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].addNeigh(&_nodelist[indexIn]); // add inNode of edge to neighbours of outNode of edge
    }
}

// function that removes an edge from the edgelist of the clustered graph + removes the corresponding neighbours of the involved nodes
void Clustered_Random_Network::removeEdge(Edge e){
    if (contains(_edgelist, e)){
        _edgelist.erase(remove(_edgelist.begin(), _edgelist.end(), e), _edgelist.end());
        int indexIn = e.inNode()->index();
        int indexOut = e.outNode()->index();
        _nodelist[indexIn].removeNeigh(&_nodelist[indexOut]); // remove outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].removeNeigh(&_nodelist[indexIn]); // remove inNode of edge to neighbours of outNode of edge
    }
}

// function that makes an Erdos-Renyi graph --> will be one of the clustered components of the clustered graph
void Clustered_Random_Network::makeErdosRenyi(int numberOfNodes, double edgeProb, int indexStart){
    random_device rd; // will be used to obtain a seed for the random number engine
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
    }
}

// function that removes all the edges from the clustered graph
void Clustered_Random_Network::removeAllEdges(){
    // remove all the corresponding neighbours
    for (int i = 0; i < _edgelist.size(); i++){
        int indexIn = _edgelist[i].inNode()->index();
        int indexOut = _edgelist[i].outNode()->index();
        _nodelist[indexIn].removeAllNeigh(); // remove neighbours of inNode of edge
        _nodelist[indexOut].removeAllNeigh(); // remove neighbours of outNode of edge
    }
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
            //Node* n = &out.back();
            //cout << n->index() << endl;
            if (&out.back() != _edgelist[i].inNode()){
                //Edge e = Edge(_edgelist[i]->inNode(), &out.back());
                //cout << e << endl;
                helpEdgelist.push_back(Edge(_edgelist[i].inNode(), &out.back())); // add new edge to the help vector
            }
            //cout << e.outNode()->index() << endl;
            //cout << helpEdgelist.back();
            //removeEdge(_edgelist[i]); // remove the old edge from the edgelist
            //--i;
        }
        else{
            helpEdgelist.push_back(_edgelist[i]);
        }
    }
    removeAllEdges();
    
    // all new made edges seem fine, except the first (always has zero as outNode...)
    for (int i = 0; i < helpEdgelist.size(); i++){
        cout << helpEdgelist[i];
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
}