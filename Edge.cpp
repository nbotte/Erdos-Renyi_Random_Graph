// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

#define _USE_MATH_DEFINES
#include "Edge.h"
#include "Node.h"
#include <iomanip>
#include <iostream>
using namespace std;

// default constructor
Edge::Edge(){};

// implementation of constructor, construct edge by defining its inNode and its outNode
Edge::Edge(shared_ptr<Node> inNode, shared_ptr<Node> outNode){
    _inNode = inNode;
    _outNode = outNode;
}

// implementation of getters, provides access to data member with corresponding name
shared_ptr<Node> Edge::inNode() const {return _inNode;}
shared_ptr<Node> Edge::outNode() const {return _outNode;}

// function that overwrites the << operator to print edges to screen
ostream& operator<<(ostream& os, const Edge& e){
    os << e._inNode->index() << '-' << e._outNode->index() << ' ';
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
