// Nina Botte

#define _USE_MATH_DEFINES
#include <iomanip>
#include <iostream>
using namespace std;

#ifndef EDGE_H
#define EDGE_H

class Node;

class Edge{  
    Node* _inNode; // declare inNode variable (pointer to startnode of edge)
    Node* _outNode; // declare outNode variable (pointer to endnode of edge)

public:
    // default constructor
    Edge();
    
    // define constructor
    Edge(Node*, Node*);

    // define copy constructor
   // Edge(const Edge &e);

    // define a destructor
 //   ~Edge();

    // define getters, provides access to data member with corresponding name
    Node* inNode() const;
    Node* outNode() const;

    // operator assignement
  //  Edge& operator= (const Edge &e);

    // operator overload << (make it a friend of Edge)
    friend ostream& operator<<(ostream& os, const Edge& e);

    // operator overload == (make it a friend of Edge)
    friend bool operator==(Edge e1, Edge e2);
};

#endif