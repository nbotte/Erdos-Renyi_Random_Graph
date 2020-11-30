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
    int index = n.index();
    _nodelist[index] = n;
}

// function that adds an edge to the edgelist of the graph + adds the corresponding neighbours to neighbourlist of the involved nodes
void Graph::addEdge(Edge e){
    // check if edge is already there
    if (contains(_edgelist, e) == false){
        //_edgelist.push_back(e); //--> takes long time?
        int indexIn = e.inNode()->index();
        int indexOut = e.outNode()->index();
        _nodelist[indexIn].addNeigh(indexOut); // add outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].addNeigh(indexIn); // add inNode of edge to neighbours of outNode of edge
    }
}

// function that removes an edge from the edgelist of the graph + removes the corresponding neighbours of the involved nodes --> NOT USED, NEEDED?
void Graph::removeEdge(Edge e){
    if (contains(_edgelist, e)){
        _edgelist.erase(remove(_edgelist.begin(), _edgelist.end(), e), _edgelist.end());
        int indexIn = e.inNode()->index();
        int indexOut = e.outNode()->index();
        _nodelist[indexIn].removeNeigh(indexOut); // remove outNode of edge to neighbours of inNode of edge
        _nodelist[indexOut].removeNeigh(indexIn); // remove inNode of edge to neighbours of outNode of edge
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
}

// function to change the opinions of the nodes in graph based on majority model
void Graph::changeOpinions(){ 
    // first: each active node sends it current opinion to its neighbours
    for (int i = 0; i < _nodelist.size(); i++){
        if (_nodelist[i].active()){
            // active, so sends its current opinion to the opinionlist of its neighbours (all of them)
            // But, the opinion send at the current timestep should not play a role in updating the opinion of active neighbours in this timestep
            // --> so set the current opinion equal to old opinion, then update the opinions and then send the old opinion 
            _nodelist[i].setOldOpinion(_nodelist[i].opinion());
        }
    }

    // second: each active node gets a new opinion according to the opinions in its opinionlist
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].changeOpinion(); // changeOpinion checks if the node is active, so only active nodes get updated
    }

    // now send the oldOpinion of the active nodes to all their neighbours
    for (int i = 0; i < _nodelist.size(); i++){
        if (_nodelist[i].active()){
            for (int index : _nodelist[i].neigh()){
                _nodelist[index].addOpinion(_nodelist[i].oldOpinion());
            }
        }
    }



    /*for (int i = 0; i < _nodelist.size(); i++){
        for (int index : _nodelist[i].neigh()){
            _nodelist[i].addNeighOpinion(_nodelist[index].opinion());
        }
    }
    // second: each active node gets a new opinion according to the opinions in its opinionlist
    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].changeOpinion();
    } */
    /*for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].setNewOpinion();
        _nodelist[i].removeAllNeighOpinion();
        if (_nodelist[i].active()){
            for (int index : _nodelist[i].neigh()){
                _nodelist[index].addOpinion(_nodelist[i].oldOpinion());
            } 
        }
    }*/
/*for (int i = 0; i < _nodelist.size(); i++){
        if (_nodelist[i].active()){
            _nodelist[i].setOldOpinion(_nodelist[i].opinion());
        }
    }*/
}

void Graph::changeRandomOpinion(){
    size_t nelem = 1;
    vector<Node> out;

    sample(_nodelist.begin(), _nodelist.end(), back_inserter(out), nelem, mt19937{random_device{}()});
    //cout << _nodelist[out.back().index()] << ": ";
    for (int index : _nodelist[out.back().index()].neigh()){
        // if the neighbour is active: add its opinion to the neighOpinion list of the current Node
        if (_nodelist[index].active()){
            _nodelist[out.back().index()].addNeighOpinion(_nodelist[index].opinion());
           // cout << _nodelist[index] << ' ';
        }
    }
   // cout << endl;
    _nodelist[out.back().index()].changeOpinion(); // changeOpinion checks if the node is active, so only active nodes can update their opinion
  //  _nodelist[out.back().index()].setNewOpinion();
    _nodelist[out.back().index()].removeAllNeighOpinion();
}

// function to deactivate all the nodes in the network, but first the current active nodes need to set their wasActive variable to true
// this function might be unneccessary!
void Graph::deactivateNodes(){
    for (int i = 0; i < _nodelist.size(); i++){
        //_nodelist[i].setWasActive();
        _nodelist[i].deactivate();
    }
}

// function that sets a fraction of the nodes as active (according to a bernouilli distribution)
void Graph::setNodesActive(double bernProb){
    double bernouilliProbability = bernProb;
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

// function to check if an edge is in vector 
bool Graph::contains(const vector<Edge> vec, Edge e){
    return find(vec.begin(), vec.end(), e) != vec.end();
}

