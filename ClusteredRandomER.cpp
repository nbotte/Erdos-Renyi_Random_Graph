// Nina Botte

#include <cmath>
#include "ClusteredRandomER.h"
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

// TO DO: implement header file + do some basic tests
Clustered_Random_Network::Clustered_Random_Network(double rewireProbability){
    _rewireProbability = rewireProbability;

    // start with this simple case
    Erdos_Renyi_Network g1 = Erdos_Renyi_Network(25, 0.01, 1., 0);
    Erdos_Renyi_Network g2 = Erdos_Renyi_Network(25, 0.1 ,1., 25);
    Erdos_Renyi_Network g3 = Erdos_Renyi_Network(25, 0.1 ,1., 50);
    Erdos_Renyi_Network g4 = Erdos_Renyi_Network(25, 0.01 ,1., 75);

    g1.print();
    g2.print();
    g3.print();
    g4.print();
    
    // reserve enough space for the vectors
    _nodelist.reserve(g1.nodelist().size() + g2.nodelist().size() + g3.nodelist().size() + g4.nodelist().size());
    _edgelist.reserve(g1.edgelist().size() + g2.edgelist().size() + g3.edgelist().size() + g4.edgelist().size());

    // add the nodes of the different ER graphs to the nodelist of the clustered graph
    for (int i = 0; i < g1.nodelist().size(); i++){
        addNode(g1.nodelist()[i]);
    }
    for (int i = 0; i < g2.nodelist().size(); i++){
        addNode(g2.nodelist()[i]);
    }
    for (int i = 0; i < g3.nodelist().size(); i++){
        addNode(g3.nodelist()[i]);
    }
    for (int i = 0; i < g4.nodelist().size(); i++){
        addNode(g4.nodelist()[i]);
    }
    // add the edges of the different ER graphs to the edgelist of the clustered graph
    for (int i = 0; i < g1.edgelist().size(); i++){
        _edgelist.push_back(g1.edgelist()[i]);
    }
    for (int i = 0; i < g2.edgelist().size(); i++){
        _edgelist.push_back(g2.edgelist()[i]);
    }
    for (int i = 0; i < g3.edgelist().size(); i++){
        _edgelist.push_back(g3.edgelist()[i]);    
    }
    for (int i = 0; i < g4.edgelist().size(); i++){
        _edgelist.push_back(g4.edgelist()[i]);
    }
    // rewire the edges
    //rewireEdges();
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
        _nodelist[indexIn].addAdjEdge(&e); // add outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].addAdjEdge(&e); // add inNode of edge to neighbours of outNode of edge
    }
}

// function that removes an edge from the edgelist of the clustered graph + removes the corresponding neighbours of the involved nodes
void Clustered_Random_Network::removeEdge(Edge e){
    if (contains(_edgelist, e)){
        _edgelist.erase(remove(_edgelist.begin(), _edgelist.end(), e), _edgelist.end());
        int indexIn = e.inNode()->index();
        int indexOut = e.outNode()->index();
        _nodelist[indexIn].removeAdjEdge(&e); // remove outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].removeAdjEdge(&e); // remove inNode of edge to neighbours of outNode of edge
    }
}

// function that removes all the edges from the clustered graph
void Clustered_Random_Network::removeAllEdges(){
    // remove all the corresponding neighbours
    for (int i = 0; i < _edgelist.size(); i++){
        int indexIn = _edgelist[i].inNode()->index();
        int indexOut = _edgelist[i].outNode()->index();
        _nodelist[indexIn].removeAllAdjEdges(); // remove neighbours of inNode of edge
        _nodelist[indexOut].removeAllAdjEdges(); // remove neighbours of outNode of edge
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
                //Edge e = Edge(_edgelist[i].inNode(), &out.back());
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

/* WEIRD SITUATION: If I print the edges of the clustered graph, I always get a pointer to the inNode if the inNode is the first node in one of the ER graphs, however when I print this edge
in the constructor I get the node... */

/* Needs further testing, also check if neighbours change accordingly, maybe implement things nicer + think about weird situation -> how to solve (also happens to other nodes after rewireing) */

// TO DO: solve issue with indexes! Leads to segmentation faults...