// use github, not github.ugent

#define _USE_MATH_DEFINES
#include <cmath>
//#include <math.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
using namespace std;

class Node;

class Edge{
    Node* _inNode; // declare inNode variable
    Node* _outNode; // declare outNode variable

public:
    // constructor, construct edge by defining its inNode and its outNode
    Edge(Node* inNode, Node* outNode){_inNode=inNode; _outNode=outNode;}

    // getter, provides access to data member with corresponding name
    Node* inNode() {return _inNode;}
    Node* outNode() {return _outNode;}

    // THIS FUNCTION DOESN'T WORK
    // print function
    /*void print(){
        cout << to_string((int) *_inNode) << '-' << to_string((int) *_outNode) << endl;
    }*/
};

class Node{
    int _index; // declare index variable (= name of node)
    vector<Edge*> _outEdges; // declare outEdges variable (= vector of edges attached to the node)
    
public:
    // constructor, construct a node by defining its name and a vector of outgoing edges
    Node(int index){_index=index; _outEdges=vector<Edge*>();}

    // getter, provides access to data member with corresponding name
    int index() {return _index;}
    vector<Edge*> outEdges() {return _outEdges;}

    // function to add an edge to the set of outgoing edges
    void addEdge(Edge* e){
        _outEdges.push_back(e);
    }

    // function to remove an edge from the set of outgoing edges
    void removeEdge(Edge* e){
        _outEdges.erase(remove(_outEdges.begin(), _outEdges.end(), e), _outEdges.end());
    }

    // function that gives the neighbours of the node
    vector<Node> neighbours(){
        vector<Node> neighbours;
        // loop over vector of outgoing edges of node
        for (int i = 0; i < _outEdges.size(); ++i){
            Edge* e = _outEdges[i];
            if (e->inNode()->_index != _index){
                neighbours.push_back(*e->inNode());
            }
            else {
                neighbours.push_back(*e->outNode());
            }
        }
        return neighbours;
    }  
};

/*class Erdos_Renyi_Network{
    int _numberOfNodes; // total number of nodes in the graph
    double _edgeProbability; // the probablility of having an edge between any pair of nodes
    vector<Node> _nodelist; // list of nodes in graph
    vector<Edge> _edgelist; // list of edges in graph

public:
    // constructor, construct graph by making the nodes and the edges with a given probability
    Erdos_Renyi_Network(int numberOfNodes, double edgeProbability){
        _numberOfNodes = numberOfNodes;
        _edgeProbability = edgeProbability;
        
        // add nodes to the graph
        for (int i = 0; i <= _numberOfNodes; i++){
            addNode(i);
        }
        // add edges to graph with a certain probability
        random_device rd; // will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
        uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < _numberOfNodes; i++){
            for (int j = i; j < _numberOfNodes; j++){
                double r = dis(gen);
                if (r < _edgeProbability){
                    Node inNode = Node(i);
                    Node outNode = Node(j);
                    addEdge(inNode, outNode);
                }
            }
        }
    }

    // function to add a node to the graph
    void addNode(int index){
        Node n = Node(index);
        _nodelist.push_back(n);
    }

    // function to add an edge to the graph
    void addEdge(Node inNode, Node outNode){
        Edge e = Edge(inNode, outNode);
        // check if edge is already there
        if (contains(_edgelist, e) == false){
            _edgelist.push_back(e);
            inNode.addEdge(e); // add this new edge to the set of edges attached to inNode
            outNode.addEdge(e); // add this new edge to the set of edges attached to outNode
        }
    }

    // function to remove an edge from the graph
    void removeEdge(Edge e){
        _edgelist.erase(e); // remove edge from edgelist
        e.inNode().removeEdge(e); // remove edge from the set of edges attached to inNode
        e.outNode().removeEdge(e); // remove edge from the set of edges ettached to outNode
    }

    // function to remove a node from the graph
    void removeNode(Node n){
        vector<Edge> edges = n.outEdges();
        // remove all edges attached to node that needs to be removed
        for (int i = 0; i <= edges.size(); ++i){
            Edge e = edges.at(i);
            removeEdge(e);
        }
        _nodelist.erase(n.index()); // remove node from nodelist NOT GOOD YET!!
    }

    // function to check if an element is in vector
    bool contains(const vector<Edge> &vec, Edge e){
	    return find(vec.begin(), vec.end(), e) != vec.end();
    }

    // print function 
    void print(){
        cout << 'Nodes: ';
        for (int i = 0; i < _nodelist.size(); i++){
            cout << _nodelist.at(i) << ' '; 
        }
        cout << endl;
        cout << 'Edges: ';
        for (int i = 0; i < _edgelist.size(); i++){
            cout << _edgelist.at(i) << ' '; 
        }
        cout << endl;
    }
};*/

int main(){
    double p = 0.01;
    int N = 100;
    Node n = Node(1);
    Node m = Node(2);
    Edge e = Edge(&n, &m);
    n.addEdge(&e);
    m.addEdge(&e);
    cout << m.index() << endl;
    // print function --> seems to work
    for (int i = 0; i < m.neighbours().size(); i++){
        cout << m.neighbours()[i].index() << ' ';
    }
    cout << endl;
    
};

// maybe also include adjecency matrix
// https://stackoverflow.com/questions/4964482/how-to-create-two-classes-in-c-which-use-each-other-as-data --> explains 
// how to use 2 classes that reference each other, maybe reconstruct with header files?
// debug network class!
// node and edge class seem to be okay, but is there a nicer way?