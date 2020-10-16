// Nina Botte

#define _USE_MATH_DEFINES
#include "Edge.h"
#include "Node.h"
#include <iomanip>
#include <iostream>
using namespace std;

// implementation of constructor, construct edge by defining its inNode and its outNode
Edge::Edge(Node* inNode, Node* outNode){_inNode=inNode; _outNode=outNode;}

// implementation of getters, provides access to data member with corresponding name
Node* const Edge::inNode() {return _inNode;}
Node* const Edge::outNode() {return _outNode;}

// function that overwrites the << operator to print edges to screen
ostream& operator<<(ostream& os, Edge& e){
    os << e.inNode()->index() << '-' << e.outNode()->index() << ' ';
    return os;
}

// function that overwrites the == operator to compare 2 edges
bool operator==(Edge e1, Edge e2){
    if (e1.inNode() == e2.inNode()){
        return e1.outNode() == e2.outNode();
    }
    else if (e1.inNode() == e2.outNode()){
        return e1.outNode() == e2.inNode();
    }
    else{ return false;}
}
