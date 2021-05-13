// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

#define _USE_MATH_DEFINES
#include <iomanip>
#include <iostream>
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
    int _oldOpinion; // declare oldOponion variable (will contain the old opinion of the node when updating the opinions)
    double _resistance; // declare resistance variable (= stubborness of the node, resistance to change his opinion)
    double _threshold; // variable that contains the treshold for a node to change opinion or not (kind of determines its stubborness)
    bool _active; // declare active variable (= determines if node is active or not)
    vector<int> _neighOpinion; // this is the hidden list (use vector dataset) with the opinions that the neighbors posted since the last time the node was active (see paper 8)
    int _cluster; // variable that determines to what cluster the node belongs; mainly used in the SBM, BUT can also be used to optimize modularity etc... 

public: 
    Node(); // define default constructor
    Node(int, int, double, bool); // define constructor

    // define getters, provides access to data member with corresponding name
    int index() const;
    list<int> neigh() const;
    list<int> helpNeigh() const;
    int opinion() const;
    int oldOpinion() const;
    double resistance() const;
    double threshold() const;
    bool active() const;
    int cluster() const;
    vector<int> neighOpinion() const;

    // declare member functions of class Node
    // functions for adding and removing neighbors
    void addNeigh(int index);
    void addHelpNeigh(int index);
    void removeNeigh(int index);
    void removeAllNeigh();
    void removeAllHelpNeigh();
    // functions for changing opinion
    void changeOpinion();
    void setOpinion(int opinion);
    void addOpinion(int opinion);
    void setOldOpinion(int opinion);
    void addNeighOpinion(int opinion);
    void removeAllNeighOpinion();
    // functions regarding being active or not
    void deactivate();
    void setActive(bool active);
    // functions regarding stubbornness
    void setResistance(double resistance);
    void setThreshold(double threshold);
    // function about cluster to which node belongs
    void setCluster(int cluster);
    // functions for ordering the opinions in the hidden time-line
    void orderOpinionsPR(); 
    void orderOpinionsREC();
    void orderOpinionsREF();
    // function to see if certain node is neighbor of current node
    bool containsNeigh(int n); 

    // operator overload == (make it a friend of Node class)
    friend bool operator==(Node n1, Node n2);

    // operator overload << (make it a friend of Node class)
    friend ostream& operator<<(ostream& os, const Node& n);
};



#endif