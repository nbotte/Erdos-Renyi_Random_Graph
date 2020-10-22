// Nina Botte

#define _USE_MATH_DEFINES
#include "Edge.h"
#include "Node.h"
#include <iomanip>
#include <iostream>
using namespace std;

// implementation of constructor, construct edge by defining its inNode and its outNode
Edge::Edge(Node* inNode, Node* outNode){
    // allocate the memory to hold a node
    _inNode = new Node;
    _outNode = new Node;

    *_inNode = *inNode; 
    *_outNode = *outNode;
}

// implementation of the copy constructor
Edge::Edge(const Edge &e){
    // allocate memory first
    _inNode = new Node;
    _outNode = new Node;

    // then copy the value from the passed object
    *_inNode = *e._inNode; 
    *_outNode = *e._outNode;
}

// destructor
Edge::~Edge(){
    // release the memory allocated
    delete _inNode;
    delete _outNode;
    _inNode = NULL;
    _outNode = NULL;
}

// implementation of getters, provides access to data member with corresponding name
Node* Edge::inNode() const {return _inNode;}
Node* Edge::outNode() const {return _outNode;}

// implementation of operator assignement
Edge & Edge::operator=(const Edge &e){
    // check for self assignement
    if (this != &e){
        *_inNode = *(e._inNode);
        *_outNode = *(e._outNode);
    }
    return *this;
}

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
