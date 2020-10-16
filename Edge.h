// Nina Botte

#define _USE_MATH_DEFINES
#include "Node.h"
#include <iomanip>
#include <iostream>
using namespace std;

#ifndef EDGE_H
#define EDGE_H

class Edge{  
    Node* _inNode; // declare inNode variable (pointer to startnode of edge)
    Node* _outNode; // declare outNode variable (pointer to endnode of edge)

public:
    // define constructor
    Edge(Node*, Node*);

    // define getters, provides access to data member with corresponding name
    Node* const inNode();
    Node* const outNode();

    // operator overload << (make it a friend of Edge)
    friend ostream& operator<<(ostream& os, Edge& e);

    // operator overload == (make it a friend of Edge)
    friend bool operator==(Edge e1, Edge e2);
};

#endif