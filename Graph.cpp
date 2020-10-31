// Nina Botte

#include <cmath>
#include <memory>
#include "Graph.h"
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

// virtual destructor
Graph::~Graph(){};

// implementation of the getters, provides access to data member with corresponding name
vector<Node> Graph::nodelist() {return _nodelist;}
vector<Edge> Graph::edgelist() {return _edgelist;} // --> Not needed anymore?

// virtual function makeGraph
void Graph::makeGraph(){};

// function to add a node to the graph
void Graph::addNode(Node n){
    _nodelist.push_back(n);
}

// function that adds an edge to the edgelist of the graph + adds the corresponding neighbours to neighbourlist of the involved nodes
void Graph::addEdge(Edge e){
    // check if edge is already there
    if (contains(_edgelist, e) == false){
        //_edgelist.push_back(e); --> takes long time
        int indexIn = e.inNode()->index();
        int indexOut = e.outNode()->index();
        _nodelist[indexIn].addNeigh(make_shared<Node>(_nodelist[indexOut])); // add outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].addNeigh(make_shared<Node>(_nodelist[indexIn])); // add inNode of edge to neighbours of outNode of edge
    }
}

// function that removes an edge from the edgelist of the graph + removes the corresponding neighbours of the involved nodes --> NOT USED, NEEDED?
void Graph::removeEdge(Edge e){
    if (contains(_edgelist, e)){
        _edgelist.erase(remove(_edgelist.begin(), _edgelist.end(), e), _edgelist.end());
        int indexIn = e.inNode()->index();
        int indexOut = e.outNode()->index();
        _nodelist[indexIn].removeNeigh(make_shared<Node>(_nodelist[indexOut])); // remove outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].removeNeigh(make_shared<Node>(_nodelist[indexIn])); // remove inNode of edge to neighbours of outNode of edge
    }
}

// function that removes all the edges from the clustered graph --> NOT USED, NEEDED?
void Graph::removeAllEdges(){
    // remove all the corresponding neighbours
    for (int i = 0; i < _edgelist.size(); i++){
        int indexIn = _edgelist[i].inNode()->index();
        int indexOut = _edgelist[i].outNode()->index();
        _nodelist[indexIn].removeAllNeigh(); // remove neighbours of inNode of edge
        _nodelist[indexOut].removeAllNeigh(); // remove neighbours of outNode of edge
    }
    // Attention: normally clear() would not affect the capacity of the vector, but not always garanteed, so in case of trouble, consider this! (reserve space back after clear)
    _edgelist.clear(); // remove the edges
    _edgelist.reserve(1000);
}

// function that rewires the edges of the graph with a certain probability --> put this in clustered graph section?
// if you would throw away edgelist, than loop over nodes + loop over each neighbour and for each neighbour change it with some probability (draw random node from nodelist) --> faster run time
void Graph::rewireEdges(){
   // vector<Edge> helpEdgelist;
   // helpEdgelist.reserve(_edgelist.size());

    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    size_t nelem = 1;
    vector<Node> out;

    for (int i = 0; i < _nodelist.size(); i++){
        if (_nodelist[i].neigh().size() > 0){
            for (shared_ptr<Node> n : _nodelist[i].neigh()){
                if (n->index() > _nodelist[i].index()){
                    double r = dis(gen); // random number that will decide if edge is rewired or not
                    if (r < _rewireProbability){
                        // make sure to compile with c++17 (than sample will not give a problem)
                        sample(_nodelist.begin(), _nodelist.end(), back_inserter(out), nelem, mt19937{random_device{}()});
                        _nodelist[i].addHelpNeigh(make_shared<Node>(out.back()));
                    }
                    else{
                        _nodelist[i].addHelpNeigh(n);
                    }
                }
            }
        }    
        _nodelist[i].removeAllNeigh();        
    }    
    for (int i = 0; i < _nodelist.size(); i++){
        for (shared_ptr<Node> n : _nodelist[i].helpNeigh()){
            _nodelist[i].addNeigh(n);
            _nodelist[n->index()].addNeigh(make_shared<Node>(_nodelist[i]));
        }
        _nodelist[i].removeAllHelpNeigh();
    }




    /*for (int i = 0; i < _edgelist.size(); i++){
        double r = dis(gen); // random number that will decide if edge is rewired or not
        if (r < _rewireProbability){
            // make sure to compile with c++17 (than sample will not give a problem)
            sample(_nodelist.begin(), _nodelist.end(), back_inserter(out), nelem, mt19937{random_device{}()});
            if (make_shared<Node>(out.back()) != _edgelist[i].inNode()){
                helpEdgelist.push_back(Edge(_edgelist[i].inNode(), make_shared<Node>(out.back()))); // add new edge to the help vector
            }
        }
        else{
            helpEdgelist.push_back(_edgelist[i]);
        }
    }
    removeAllEdges();

    for (int i = 0; i < helpEdgelist.size(); i++){
        addEdge(helpEdgelist[i]); // add the new edges in help vector to the edgelist of the clustered graph
    }*/
}

// function to change the opinions of the nodes in graph based on majority model
void Graph::changeOpinions(){
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
void Graph::deactivateNodes(){
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].setWasActive();
        _nodelist[i].deactivate();
    }
}

// function that sets a fraction of the nodes as active (according to a bernouilli distribution)
void Graph::setNodesActive(){
    double bernouilliProbability = 1.;
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    bernoulli_distribution disBern(bernouilliProbability);
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].setActive(disBern(gen));
    }
}

// function that counts the fraction of the 2 opinions in the graph (returns a vector with the 2 fractions)
vector<double> Graph::countOpinionFraction(){
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
    fractions.push_back(double(opinion0)/double(_nodelist.size()));
    fractions.push_back(double(opinion1)/double(_nodelist.size()));
    return fractions;
}

// print function for a graph --> prints nodes and edges
void Graph::print(){
    cout << "Nodes: ";
    for (int i = 0; i < _nodelist.size(); ++i){
        cout << _nodelist[i]; 
    }
    cout << endl;
        
    cout << "Edges: "; // No edges anymore?
    for (int i = 0; i < _edgelist.size(); i++){
        cout << _edgelist[i]; 
    }
    cout << endl;
}

// function to check if an element is in vector 
bool Graph::contains(const vector<Edge> vec, Edge e){
    return find(vec.begin(), vec.end(), e) != vec.end();
}