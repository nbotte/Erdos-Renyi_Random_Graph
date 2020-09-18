// use github, not github.ugent

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <time.h>
using namespace std;

class Node{

    // represents the status of the node: -1 is empty, >=0 is index of node
    int _status;

    // pointer to edges attached to that node
    vector<Edge*> _EdgesOfNode; 
    // QUESTION: I construct a node, but I don't have edges yet, how can I point to the edges then?

    // pointer to neighbours of node
    vector<Node*> _Neighbours;
    // QUESTION: How to determine neighbours of a node?

public:
    // constructor of node, default node is empty (status -1) with no edges and no neighbours
    Node(){_status = -1;}

    // DESTRUCTOR NEEDED?

    // basis public functions
    int status() const {return _status;}
    vector<Edge*> EdgesOfNode() const {return _EdgesOfNode;}    
    vector<Node*> Neighbours() const {return _Neighbours;}

    // functions to add and remove node (= change status from or to empty)
    // takes as argument the index of the node you want to add or remove
    void add(int index){
        if (_status == -1){
            _status = index;
        }
    }
    void remove(int index){
        if (_status == index){
            _status = -1;
        }
    }
};

class Edge{
    // needs pointer to nodes it connect
    // needs function that returns those nodes
    // needs function that adds and removes an edge

};

class Erdos_Renyi_Network{
    // object of this class will be an edge list and a node list
    // needs function to construct N nodes
    // needs function to create edges with probability p (faster than N^2 ?)

};