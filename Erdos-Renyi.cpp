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
};

class Node{
    int _index; // declare index variable (= name of node)
    vector<int> _neigh; // declare neigh variable (= vector of indices of nodes neighbouring the current node)
    
public:
    // constructor, construct a node by defining its name and a vector of neighbours
    Node(int index){_index=index; _neigh=vector<int>();}

    // getter, provides access to data member with corresponding name, should const be there?
    int index() {return _index;}
    vector<int> neigh() {return _neigh;}

    // function to add a neighbour to the vector of neighbours of a node
    void addNeigh(int index){
        _neigh.push_back(index);
    }

    // function to remove a neighbour from the vector of neighbours of a node, NOT TESTED
    void removeNeigh(int index){
        _neigh.erase(remove(_neigh.begin(), _neigh.end(), index), _neigh.end());
    }

    // function that gives the neighbours of the node STILL NECESSARY?
    /*vector<Node*> neighbours(){
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
    vector<Node> _nodelist; // list of nodes in graph
    //vector<vector<int>> _neighlist; // list of neighbours of each node in graph
    vector<Edge> _edgelist; // list of edges in graph

public:
    // constructor, construct graph by making the nodes and the edges with a given probability
    Erdos_Renyi_Network(int numberOfNodes, double edgeProbability){
        _numberOfNodes = numberOfNodes;
        _edgeProbability = edgeProbability;

        _nodelist.reserve(_numberOfNodes); // needed?
        _edgelist.reserve(pow(_numberOfNodes, 2)); //needed?
        // add nodes to the graph 
        for (int i = 0; i < _numberOfNodes; i++){
            Node n = Node(i);
            addNode(n);
        }    
        // add edge between any pair of nodes with a certain probability
        random_device rd; // will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
        uniform_real_distribution<> dis(0.0, 1.0);
        //vector<int> neighbours;
        for (int i = 0; i < _nodelist.size(); i++){
            for (int j = i; j < _nodelist.size(); j++){
                double r = dis(gen);
                if (r < _edgeProbability){
                    addEdge(i, j);
                    //_nodelist[i].addNeigh(_nodelist[j].index()); // add this new edge to the set of edges attached to inNode
                    //_nodelist[j].addNeigh(_nodelist[i].index()); // add this new edge to the set of edges attached to outNode
                   // neighbours.push_back(_nodelist[j].index()); // maybe put in function add neighbour?
                }
            }
            //_neighlist.push_back(neighbours);
            //neighbours.clear();
        }
    }

    // getter, provides access to data member with corresponding name
    vector<Node> nodelist() {return _nodelist;}
    //vector<vector<int>> neighlist() {return _neighlist;} 
    vector<Edge> const edgelist() {return _edgelist;}

    // function to add a node to the graph
    void addNode(Node n){
        _nodelist.push_back(n);
    }

    // function to add an edge to the graph
    void addEdge(int indexIn, int indexOut){
        Node* N = new Node(indexIn);
        Node* M = new Node(indexOut);
        Edge e = Edge(N, M);
        // check if edge is already there
        if (contains(_edgelist, e) == false){
            _edgelist.push_back(e);
            _nodelist[indexIn].addNeigh(_nodelist[indexOut].index()); // add outNode of edge to neighbours of inNode of edge
            _nodelist[indexOut].addNeigh(_nodelist[indexIn].index()); // add inNode of edge to neighbours of outNode of edge

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
        vector<int> neigh = _nodelist[index].neigh();
        // remove all edges attached to node that needs to be removed
        for (int i = 0; i < neigh.size(); ++i){
            Node m = neigh[i];
            Edge* e = new Edge(n, &m);
            removeEdge(e);
        }
        _nodelist.erase(&_nodelist[index]); // remove node from nodelist NOT GOOD YET!!
    }*/

    // function that implements a comparison of 2 edges (== operator), returns a bool
    static bool equalEdge(Edge e1, Edge e2){
        if((e1.inNode()->index() == e2.inNode()->index() && e1.outNode()->index() == e2.outNode()->index()) ||
        (e1.inNode()->index() == e2.outNode()->index() && e1.outNode()->index() == e2.inNode()->index())){
            return true;
        }
        else{
            return false;
        }
    }

    // function to check if an element is in vector
    bool contains(const vector<Edge> vec, Edge e){
	    return find_if(vec.begin(), vec.end(), bind(equalEdge, std::placeholders::_1, e)) != vec.end();
    }

    // print function 
    void print(){
        cout << "Nodes: ";
        for (int i = 0; i < _nodelist.size(); ++i){
            cout << _nodelist[i].index() << ' '; 
        }
        cout << endl;
        /*cout << "Neighbours: ";
        for (int i = 0; i < _neighlist.size(); i++){
            cout << "Node: " << i << ": ";
            for (int j = 0; j < _neighlist[i].size(); j++){
                cout << _neighlist[i][j] << ' ';
            }
            cout << endl; 
        }*/
        cout << "Edges: ";
        for (int i = 0; i < _edgelist.size(); i++){
            cout << _edgelist[i].inNode()->index() << '-' << _edgelist[i].outNode()->index() << ' '; 
        }
        cout << endl;
    }
};

int main(){
    // This example seem to work as expected
   /* Node n = Node(1);
    Node m = Node(2);
    Node k = Node(0);
    Node l = Node(3);

    Edge e1 = Edge(&n, &m);
    Edge e2 = Edge(&n, &k);
    Edge e3 = Edge(&n, &l);
    Edge e4 = Edge(&m, &k);
    Edge e5 = Edge(&l, &k);

    n.addNeigh(m.index());
    m.addNeigh(n.index());
    n.addNeigh(k.index());
    k.addNeigh(n.index());
    n.addNeigh(l.index());
    l.addNeigh(n.index());
    m.addNeigh(k.index());
    k.addNeigh(m.index());
    l.addNeigh(k.index());
    k.addNeigh(l.index());

    n.removeNeigh(m.index());

    vector<int> neigh = n.neigh();
    for (int i = 0; i < neigh.size(); i++){
        cout << neigh[i] << endl;   
    }*/

    // SEEMS TO BE FINE!! :o
    int N = 100;
    double p = 0.01;
    Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p);
    g.print();
    // WORKS!
    for (int i = 0; i < g.nodelist()[0].neigh().size(); i++){
        cout << g.nodelist()[0].neigh()[i] << ' ';
    }
};

// maybe also include adjecency matrix
// https://stackoverflow.com/questions/4964482/how-to-create-two-classes-in-c-which-use-each-other-as-data --> explains 
// how to use 2 classes that reference each other, maybe reconstruct with header files?
// Not a mess anymore, declaring new pointers by 'new' seems to be the trick! However, they should also be manually be removed from memory, where in the code should this be done? In destructor?

// Problem: if a make a graph, the function addNeigh of Node class doesn't behave properly, problem with use of pointers? Or problem in constructor of node?
// Possible solution: remove edge class, only work with nodes and their neighbours, avoid use of pointers?

// current status: edge class is used anymore, but is it necessary/valuable? Adding neighbours seems to work fine! Remove node and remove edge from graph are not ok yet, but not sure if they are necessary
// TO DO: node should have some extra properties: active, stubborn, opinion
// TO DO: opinion dynamics over time, add method to graph class to change opinions of nodes based on some condition