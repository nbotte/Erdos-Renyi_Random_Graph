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
    // this class should also contain a function initialize_opinion that initializes the opinion (A or B) of the particular node (randomly)
    // also a function opinion is needed that determines the opinion of the node based on the opinion of the neighbours
    // Question: should these functions be incorporated in this class or should I make a different class that deals with the opinion dynamics?

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

    // status = false: no edge, status = true: edge
    bool _status;
    // start and end node of an edge
    Node _start;
    Node _end;

    // pointer to the pair of nodes that the edge connects
    vector<Node *> _PairOfNodes;

public:
    // constructor of the edge (default edge is non-existing, status = false)
    Edge(){_status = false}    

    // basic public functions
    bool status() const {return _status;}
    Node start() const {return _start;}
    Node end() const {return _end;}
    vector<Node *> PairOfNodes() const {return _PairOfNodes;}

    // functions that adds/removes an edge between two nodes
    void add(Node v, Node w){
        if (_status == false){
            _status = true;
            _start = v;
            _end = w;
            _PairOfNodes.push_back(_start);
            _PairOfNodes.push_back(_end);
        }
    void remove(Node v, Node w){
        if (_status == true){
            _status = false;
            _PairOfNodes.clear();
        }
    }

};

class Erdos_Renyi_Network{
    // object of this class will be an edge list and a node list
    // needs function to construct N nodes
    // needs function to create edges with probability p (faster than N^2 ?)

    vector<Node> _NodeList;
    vector<Edge> _EdgeList;
    int _NumberOfNodes;
    double _EdgeProbability;

public: 
    // constructor of graph
    Erdos_Renyi_Network(const double EdgeProbability, const int NumberOfNodes){
        _EdgeProbability = EdgeProbability;
        _NumberOfNodes = NumberOfNodes;
        for (int i = 0; i <= NumberOfNodes - 1; i++){
            Node v;
            v.add(i);
            _NodeList.push_back(v);
        }
        // for now all the edges are added, just as as test!
        for (int n = 0; n <= NumberOfNodes - 1; n++){
            for (int m = 0; m <= NumberOfNodes - 1; m++){
                if (n != m){
                    Edge e;
                    e.add(n, m);
                    _EdgeList.push_back(e);
                }
            }
        }
    }

    // basic public functions
    vector<Node> NodeList() const {return _NodeList}
    vector<edge> EdgeList() const {return _EdgeList}

};

int main(){
    double p = 0.1;
    int N = 100;
    Erdos_Renyi_Network graph = Erdos_Renyi_Network(p, N);
    edges = graph.EdgeList();
    nodes = graph.NodesList();

    cout << edges;
    cout << nodes;

}