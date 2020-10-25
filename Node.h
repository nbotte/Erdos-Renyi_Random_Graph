// Nina Botte

#define _USE_MATH_DEFINES
#include <iomanip>
#include <iostream>
#include <list>
using namespace std;

#ifndef NODE_H
#define NODE_H

class Node{  
    int _index; // declare index variable (= name of node)
    list<Node*> _neigh; // declare neigh variable (= vector of pointers to nodes neighbouring the current node)
    int _opinion; // declare opinion variable (= opinion of node at current time step, choice between 0 and 1)
    int _newOpinion; // declare newOpinion variable (= opinion of node at the next time step, choice between 0 and 1)
    double _resistance; // declare resistance variable (= stubborness of the node, resistance to change his opinion)
    bool _active; // declare active variable (= determines if node is active or not)
    bool _wasActive; // declare wasActive variable (= determines if the node was active in the previous timestep)

public: 
    Node(); // define default constructor
    Node(int, int, double, bool); // define constructor

    // define getters, provides access to data member with corresponding name
    int index() const;
    list<Node*> neigh() const;
    int opinion() const;
    int newOpinion() const;
    double resistance() const;
    bool active() const;

    // declare member functions of class Node
    void addNeigh(Node* n);
    void removeNeigh(Node* n);
    void removeAllNeigh();
    void changeOpinion();
    void setNewOpinion();
    void deactivate();
    void setActive(bool active);
    void setWasActive();

    // operator overload == (make it a friend of Node)
    friend bool operator==(Node n1, Node n2);

    // operator overload << (make it a friend of Node)
    friend ostream& operator<<(ostream& os, const Node& n);
};



#endif