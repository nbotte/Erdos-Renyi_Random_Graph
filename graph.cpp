#include <iostream>
#include <random>
#include <list>
#include <algorithm>
#include <set>
using namespace std;

class Node;

class Edge{
    Node* _inNode; // declare inNode variable
    Node* _outNode; // declare outNode variable

public:
    // constructor, construct edge by defining its inNode and its outNode
    Edge(Node* inNode, Node* outNode){_inNode=inNode; _outNode=outNode;}

    // getter, provides access to data member with corresponding name
    Node* const inNode() {return _inNode;}
    Node* const outNode() {return _outNode;}
};

bool operator==(Edge e1, Edge e2){
    if (e1.inNode() == e2.inNode()){
        return e1.outNode() == e2.outNode();
    }
    else if (e1.inNode() == e2.outNode()){
        return e1.outNode() == e2.inNode();
    }
    else{ return false;}
}


class  Node{
    int _name;
    set<Edge> _adjEdgelist;

public:
    Node(int name){_name=name; _adjEdgelist={};}

    int name() {return _name;}
    set<Edge> adjEdgelist() {return _adjEdgelist;}

    void addEdge(Edge e){
        _adjEdgelist.insert(e);
    }

    void removeEdge(Edge e){
        _adjEdgelist.erase(remove(_adjEdgelist.begin(), _adjEdgelist.end(), e), _adjEdgelist.end());
    }
};

/*class Graph{
    int _numberOfNodes;
    double _edgeProbability;
    list<Node> _nodelist;

public:
    Graph(int numberOfNodes, double edgeProbability){
        _numberOfNodes = numberOfNodes;
        _edgeProbability = edgeProbability;

        random_device rd; // will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
        uniform_real_distribution<> dis(0.0, 1.0);

        for (int i = 0; i < _numberOfNodes; i++){
            Node n = Node(i);
            addNode(n);
        }
        for (int i = 0; i < _nodelist.size(); i++){
            for (int j = i; j < _nodelist.size(); j++){
                double r = dis(gen);
                if (r < _edgeProbability){
                    addEdge(i, j);
                }
            }
        }
    } 

    list<Node> nodelist() {return _nodelist;}

    void addNode(Node n){
        _nodelist.push_back(n);
    }   

    void addEdge(Node* n, Node* m){
        Edge e = Edge(n, m);
        // check if edge is already there
        if (contains(_edgelist, e) == false){
            _edgelist.push_back(e);
            _nodelist[indexIn].addNeigh(&_nodelist[indexOut]); // add outNode of edge to neighbours of inNode of edge
            _nodelist[indexOut].addNeigh(&_nodelist[indexIn]); // add inNode of edge to neighbours of outNode of edge

        }
    }
};*/

ostream& operator<<(ostream& os, Edge& e){
    os << e.inNode()->name() << '-' << e.outNode()->name() << ' ';
    return os;
}

int main(){
    Node n = Node(1);
    Node m = Node(2);
    Node k = Node(4);
    Edge e1 = Edge(&n, &m);
    Edge e2 = Edge(&n, &k);

    n.addEdge(e1);
    for (Edge e : n.adjEdgelist()){
        cout << e;
    }
    /*n.removeEdge(e1);
    for (const auto &e : n.adjEdgelist()){
        cout << e;
    }*/

    cout << (e1==e2) << endl;
    cout << (e1==e1) << endl;

}