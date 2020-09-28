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
#include <functional>
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
    vector<int> _neigh; // declare neigh variable (= vector of indices of nodes neighbouring the current node)
    
public:
    // constructor, construct a node by defining its name and a vector of outgoing edges
    Node(int index){_index=index; _neigh=vector<int>();}

    // getter, provides access to data member with corresponding name
    int index() {return _index;}
    vector<int> neigh() {return _neigh;}

    // function to add an edge to the set of outgoing edges
    void addNeigh(int index){
        _neigh.push_back(index);
        cout << _neigh.size() << ' '; // always gives size equal to 1
    }

    // function to remove an edge from the set of outgoing edges NEEDS TO BE CHANGED
    /*void removeOutEdge(Edge* e){
        _outEdges.erase(remove(_outEdges.begin(), _outEdges.end(), e), _outEdges.end());
    }

    // function that gives the neighbours of the node STILL NECESSARY?
    vector<Node*> neighbours(){
        vector<Node*> neighbours;
        // loop over vector of outgoing edges of node
        for (int i = 0; i < _outEdges.size(); ++i){
            if (_outEdges[i]->inNode()->_index != _index){
                neighbours.push_back(_outEdges[i]->inNode());
            }
            else {
                neighbours.push_back(_outEdges[i]->outNode());
            }
        }
        return neighbours;
    }  */
};


class Erdos_Renyi_Network{
    int _numberOfNodes; // total number of nodes in the graph
    double _edgeProbability; // the probablility of having an edge between any pair of nodes
    vector<Node*> _nodelist; // list of nodes in graph
    vector<Edge*> _edgelist; // list of edges in graph

public:
    // constructor, construct graph by making the nodes and the edges with a given probability
    Erdos_Renyi_Network(int numberOfNodes, double edgeProbability){
        _numberOfNodes = numberOfNodes;
        _edgeProbability = edgeProbability;

        _nodelist.reserve(_numberOfNodes); // needed?
        _edgelist.reserve(pow(_numberOfNodes, 2)); //needed?
        // add nodes to the graph 
        for (int i = 0; i < _numberOfNodes; i++){
            Node* n = new Node(i);
            addNode(n);
        }    
        // add edges to graph with a certain probability
        random_device rd; // will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
        uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < _numberOfNodes; i++){
            for (int j = i; j < _numberOfNodes; j++){
                double r = dis(gen);
                Node* inNode = new Node(i);
                Node* outNode = new Node(j);
                if (r < _edgeProbability){
                    addEdge(inNode, outNode);
                }
            }
        }
    }

    // getter, provides access to data member with corresponding name
    vector<Node*> const nodelist() {return _nodelist;}
    vector<Edge*> const edgelist() {return _edgelist;}

    // function to add a node to the graph
    void addNode(Node* n){
        _nodelist.push_back(n);
        //cout << _nodelist[n->index()]->index() << ' ';
    }

    // function to add an edge to the graph
    void addEdge(Node* n, Node* m){
        Edge* e = new Edge(n, m);
        // check if edge is already there
        if (contains(_edgelist, e) == false){
            _edgelist.push_back(e);
            Node N = *n;
            Node M = *m;
            N.addNeigh(N.index()); // add this new edge to the set of edges attached to inNode
            M.addNeigh(M.index()); // add this new edge to the set of edges attached to outNode
        }
    }

    // function to remove an edge from the graph
    // NOT TESTED YET
   /* void removeEdge(Edge* e){
        _edgelist.erase(remove(_edgelist.begin(), _edgelist.end(), e), _edgelist.end()); // remove edge from edgelist
        e->inNode()->removeOutEdge(e); // remove edge from the set of edges attached to inNode
        e->outNode()->removeOutEdge(e); // remove edge from the set of edges ettached to outNode
    }

    // function to remove a node from the graph
    // NOT TESTED YET
    void removeNode(Node* n){
        vector<Edge*> edges = n->outEdges();
        // remove all edges attached to node that needs to be removed
        for (int i = 0; i <= edges.size(); ++i){
            Edge* e = edges[i];
            removeEdge(e);
        }
        _nodelist.erase(remove(_nodelist.begin(), _nodelist.end(), n), _nodelist.end()); // remove node from nodelist NOT GOOD YET!!
    }*/

    // function that implements a comparison of 2 edges (== operator), returns a bool
    static bool equalEdge(Edge* e1, Edge* e2){
        if((e1->inNode()->index() == e2->inNode()->index() && e1->outNode()->index() == e2->outNode()->index()) ||
        (e1->inNode()->index() == e2->outNode()->index() && e1->outNode()->index() == e2->inNode()->index())){
            return true;
        }
        else{
            return false;
        }
    }

    // function to check if an element is in vector
    bool contains(const vector<Edge*> vec, Edge* e){
	    return find_if(vec.begin(), vec.end(), bind(equalEdge, std::placeholders::_1, e)) != vec.end();
    }

    // print function 
    void print(){
        cout << "Nodes: ";
        for (int i = 0; i < _nodelist.size(); ++i){
            cout << _nodelist[i]->index() << ' '; 
        }
        cout << endl;
        cout << "Edges: ";
        for (int i = 0; i < _edgelist.size(); i++){
            cout << _edgelist[i]->inNode()->index() << '-' << _edgelist[i]->outNode()->index() << ' '; 
        }
        cout << endl;
    }
};

int main(){
    // This example seem to work as expected
    /*Node n = Node(1);
    Node m = Node(2);
    Node k = Node(0);
    Node l = Node(3);

    Edge e1 = Edge(&n, &m);
    Edge e2 = Edge(&n, &k);
    Edge e3 = Edge(&n, &l);
    Edge e4 = Edge(&m, &k);
    Edge e5 = Edge(&l, &k);

    n.addOutEdge(&e1);
    m.addOutEdge(&e1);
    n.addOutEdge(&e2);
    k.addOutEdge(&e2);
    n.addOutEdge(&e3);
    l.addOutEdge(&e3);
    m.addOutEdge(&e4);
    k.addOutEdge(&e4);
    l.addOutEdge(&e5);
    k.addOutEdge(&e5);

    vector<Node> neigh = n.neighbours();
    for (int i = 0; i < neigh.size(); i++){
        cout << neigh[i].index() << endl;   
    }*/

    // SEEMS TO BE FINE!! :o
    int N = 100;
    double p = 0.01;
    Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p);
    //g.print();

    // this part is not okay, problem with neighbours
    vector<Node*> nodes = g.nodelist();
    Node n = *nodes[3];
    cout << n.neigh().size() << ' '; // gives always size equal to 0, so there must be a problem with the construction with vector of neighbours when I build the graph
    /*for (int i = 0; i < nodes.size(); i++){
        cout << nodes[0]->neigh()[i] << ' ';
    }*/
    //cout << nodes[0]->neighbours()[0]->index();
};

// maybe also include adjecency matrix
// https://stackoverflow.com/questions/4964482/how-to-create-two-classes-in-c-which-use-each-other-as-data --> explains 
// how to use 2 classes that reference each other, maybe reconstruct with header files?
// Not a mess anymore, declaring new pointers by 'new' seems to be the trick! However, they should also be manually be removed from memory, where in the code should this be done? In destructor?

// TO DO: test neighbours when you make nodes and edges individually + make/test neighbour remove function in node class + test remove node/edge functions in graph class