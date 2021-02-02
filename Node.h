// Nina Botte

#define _USE_MATH_DEFINES
#include <iomanip>
#include <iostream>
#include "boost/circular_buffer.hpp"
#include <list>
#include <vector>
using namespace std;

#ifndef NODE_H
#define NODE_H

class Node{  
    int _index; // declare index variable (= name of node)
    list<int> _neigh; // declare neigh variable (= list of indices (names) of nodes that are neighbours of the current node)
    list<int> _helpNeigh; // variable helpNeigh (will be used to rewire edges for clustered graphs)
    int _opinion; // declare opinion variable (= opinion of node at current time step, choice between 0 and 1)
    int _oldOpinion;
    int _newOpinion; // declare newOpinion variable (= opinion of node at the next time step, choice between 0 and 1)
    boost::circular_buffer<int> _opinionlist; // this list contains up to the 20 newest previous opinions of the neighbours of the node
    double _resistance; // declare resistance variable (= stubborness of the node, resistance to change his opinion)
    bool _active; // declare active variable (= determines if node is active or not)
    //bool _wasActive; // declare wasActive variable (= determines if the node was active in the previous timestep)

    vector<int> _neighOpinion; // this is the hidden list (use vector dataset) with the opinions that the neighbors posted since the last time the node was active (see paper 8)

public: 
    Node(); // define default constructor
    Node(int, int, double, bool); // define constructor

    // define getters, provides access to data member with corresponding name
    int index() const;
    list<int> neigh() const;
    list<int> helpNeigh() const;
    int opinion() const;
    int oldOpinion() const;
    int newOpinion() const;
    boost::circular_buffer<int> opinionlist() const;
    double resistance() const;
    bool active() const;

    vector<int> neighOpinion() const;

    // declare member functions of class Node
    void addNeigh(int index);
    void addHelpNeigh(int index);
    void removeNeigh(int index);
    void removeAllNeigh();
    void removeAllHelpNeigh();
    void changeOpinion();
    void setOpinion(int opinion);
    void addOpinion(int opinion);
    void setOldOpinion(int opinion);
    //void sendOpinion();
    void deactivate();
    void setActive(bool active);
    void setWasActive();
    void setResistance(double resistance);

    void addNeighOpinion(int opinion);
    void removeAllNeighOpinion();
    void orderOpinionsPR(); 
    void orderOpinionsREC();

    bool containsNeigh(int n); 

    // operator overload == (make it a friend of Node)
    friend bool operator==(Node n1, Node n2);

    // operator overload << (make it a friend of Node)
    friend ostream& operator<<(ostream& os, const Node& n);
};



#endif