// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

#define _USE_MATH_DEFINES
#include <iomanip>
#include <memory>
#include <iostream>
using namespace std;

#ifndef EDGE_H
#define EDGE_H

class Node;

class Edge{  
    shared_ptr<Node> _inNode; // declare inNode variable (pointer to startnode of edge)
    shared_ptr<Node> _outNode; // declare outNode variable (pointer to endnode of edge)

public:
    // default constructor
    Edge();
    
    // define constructor
    Edge(shared_ptr<Node>, shared_ptr<Node>);

    // define getters, provides access to data member with corresponding name
    shared_ptr<Node> inNode() const;
    shared_ptr<Node> outNode() const;

    // operator overload << (make it a friend of Edge class)
    friend ostream& operator<<(ostream& os, const Edge& e);

    // operator overload == (make it a friend of Edge class)
    friend bool operator==(Edge e1, Edge e2);
};

#endif