// Nina Botte

#include "Node.h"
#include <cmath>
#include <math.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <random>
#include <functional>
#include <algorithm>
using namespace std;

// implementation of constructor, construct a node by defining its name, opinion, resistance and a vector of neighbours (empty at construction)
Node::Node(int index, int opinion, double resistance, bool active){_index=index; _neigh=list<Node*>(); _opinion=opinion; _newOpinion=opinion; _resistance=resistance; _active=active; _wasActive=true;}

// getters
int const Node::index() {return _index;}
list<Node*> const Node::neigh() {return _neigh;}
int const Node::opinion() {return _opinion;}
int const Node::newOpinion() {return _newOpinion;}
double const Node::resistance() {return _resistance;}
bool const Node::active() {return _active;}

// function to add a neighbour to the vector of neighbours of a node
void Node::addNeigh(Node* n){
   _neigh.push_back(n);
}

// function to remove a neighbour from the vector of neighbours of a node, NOT TESTED
void Node::removeNeigh(Node* n){
    _neigh.erase(remove(_neigh.begin(), _neigh.end(), n), _neigh.end());
}

// function to change the opinion of the node
void Node::changeOpinion(){
    int opinion0 = 0; // counter that counts the number of neighbours with opinion 0
    int opinion1 = 0; // counter that counts the number of neighbours with opinion 1

    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    double r = dis(gen); // generate a random number that will determine if the resistant node will change his opinion or not
 
    // change the opinion of the active node according to the majority model and if the random number is bigger than the resistance of the node
    if (_active){
        // count the number of opinions 0 and 1 of the neighbours of the node (only take previously active neighbours into account)
        for (Node* n : _neigh){
            if (n->_wasActive){
                if (n->_opinion == 0){
                    opinion0++; 
                }
                else{ 
                    opinion1++;
                }
            }
        }

        if (opinion0 > opinion1){
            if (r >= _resistance){
                _newOpinion = 0;
            }
            else{
                _newOpinion = _opinion;
            }
        }
        else if (opinion1 > opinion0){
            if (r >= _resistance){
                _newOpinion = 1;
            }
            else{
                _newOpinion = _opinion;
            }
        }
        else{
            _newOpinion = _opinion;
        }
    }
        
}

// function that sets the opinion of a node equal to its new opinion
void Node::setNewOpinion(){
    _opinion = _newOpinion;
}

// function that deactivates the node
void Node::deactivate(){
    _active = false;
}

// function that sets the activeness of a node
void Node::setActive(bool active){
    _active = active;
}

// function that sets the wasActiveness of a node
void Node::setWasActive(){
    _wasActive = _active;
}

// overload << operator to print Node (+ opinion and activeness) to screen
ostream& operator<<(ostream& os, Node& n){
    os << n.index() << '-' << n.opinion() << '-' << n.active() << ' ';
    return os;
}


