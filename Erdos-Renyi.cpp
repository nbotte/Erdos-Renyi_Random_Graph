// use github, not github.ugent

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <random>
#include <algorithm>
using namespace std;

class Node{
    int _index; // declare index variable (= name of node)
    unordered_set<Edge> _outEdges; // declare outEdges variable (= set of edges attached to the node)

public:
    // constructor, construct a node by defining its name and a set of outgoing edges
    Node(int index){_index=index; _outEdges=unordered_set<Edge>();}

    // getter, provides read-only access to data member with corresponding name
    int index() const {return _index;}

    // function to add an edge to the set of outgoing edges
    void addEdge(Edge e){
        _outEdges.insert(e);
    }

    // function to remove an edge from the set of outgoing edges
    void removeEdge(Edge e){
        _outEdges.erase(e);
    }

    // function that gives the neighbours of the node
    unordered_set<Node> neighbours(){
        unordered_set<Node> neighbours;
        // loop over set of outgoing edges of node
        for (auto i = _outEdges.begin(); i != _outEdges.end(); ++i){
            Edge e = *i;
            if (e._inNode.index() != _index){
                neighbours.insert(e._inNode);
            }
            else {
                neighbours.insert(e._outNode);
            }
        }
        return neighbours;
    }

    // print function
    void print(Node n){
        cout << n.index() << " - Neighbours: " << n.neighbours() << endl;
    }
};

class Edge{
    Node _inNode; // declare inNode variable
    Node _outNode; // declare outNode variable

public:
    // constructor, construct edge by defining its inNode and its outNode
    Edge(Node inNode, Node outNode){_inNode=inNode; _outNode=outNode;}

    // getter, provides read-only access to data member with corresponding name
    Node inNode() const {return _inNode;}
    Node outNode() const {return _outNode;}

    // print function
    void print(Edge e){
        cout << e.inNode() << '-' << e.outNode() << endl;
    
};

class Erdos_Renyi_Network{
    int _numberOfNodes; // total number of nodes in the graph
    double _edgeProbability; // the probablility of having an edge between any pair of nodes
    list<Node> _nodelist; // list of nodes in graph
    list<Edge> _edgelist; // list of edges in graph

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
        for (int i = 0; i <= _numberOfNodes; i++){
            for (int j = i; j <= _numberOfNodes; j++){
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
        _nodelist.insert(n);
    }

    // function to add an edge to the graph
    void addEdge(Node inNode, Node outNode){
        Edge e = Edge(inNode, outNode);
        // check if edge is already there
        if (contains(_edgelist, e) == false){
            _edgelist.insert(e);
            inNode.addEdge(e); // add this new edge to the set of edges attached to inNode
            outNode.addEdge(e); // add this new edge to the set of edges attached to outNode
        }
    }

    // function to remove an edge from the graph
    void removeEdge(Edge e){
        _edgelist.remove(e); // remove edge from edgelist
        e._inNode.removeEdge(e); // remove edge from the set of edges attached to inNode
        e._outNode.removeEdge(e); // remove edge from the set of edges ettached to outNode
    }

    // function to remove a node from the graph
    void removeNode(Node n){
        unordered_set<Edge> edges = n._outEdges;
        // remove all edges attached to node that needs to be removed
        for (auto i = edges.begin(); i != edges.end(); ++i){
            Edge e = *i;
            removeEdge(e);
        _nodelist.remove(n); // remove node from nodelist
    }

    // function to check if an element is in list
    bool contains(const list<Edge> &list, Edge e){
	    return find(list.begin(), list.end(), e) != list.end();
    }
};

int main(){
    double p = 0.1;
    int N = 100;
    Erdos_Renyi_Network graph = Erdos_Renyi_Network(p, N);
    vector<edge> edges = graph.EdgeList();
    vector<node> nodes = graph.NodeList();

    cout << edges;
    cout << nodes;

};

// maybe also include adjecency matrix