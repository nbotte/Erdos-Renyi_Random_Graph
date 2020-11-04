// Nina Botte

#include "Node.h"
#include "Edge.h"
#include <iomanip>
#include <memory>
#include <iostream>
#include "boost/circular_buffer.hpp"
#include <list>
#include <random>
#include <algorithm>
using namespace std;

// attention: try not to copy nodes, becauses no proper copy constructor is written (and I do not see how to do this with neighbours that point to nodes itself)
// how to implement copy constructor?

// default constructor
Node::Node(){};

// implementation of constructor, construct a node by defining its name, opinion, resistance and a vector of neighbours (empty at construction)
Node::Node(int index, int opinion, double resistance, bool active){
    _index = index; 
    _neigh = list<shared_ptr<Node>>(); 
    _helpNeigh = list<shared_ptr<Node>>();
    _opinion = opinion; 
   // _newOpinion = opinion; 
    _opinionlist.set_capacity(20); // only keep the 20 newest opinions of the neighbours
    _resistance = resistance; 
    _active = active; 
   // _wasActive = true;
}

// implementation of the getters
int Node::index() const {return _index;}
list<shared_ptr<Node>> Node::neigh() const {return _neigh;}
list<shared_ptr<Node>> Node::helpNeigh() const {return _helpNeigh;}
int Node::opinion() const {return _opinion;}
//int Node::newOpinion() const {return _newOpinion;}
boost::circular_buffer<int> Node::opinionlist() const {return _opinionlist;}
double Node::resistance() const {return _resistance;}
bool Node::active() const {return _active;}

// function to add a neighbour to the list of pointers to the neighbours of a node
void Node::addNeigh(shared_ptr<Node> n){
   _neigh.push_back(n);
}

// function to add a help-neighbour when rewiring the clustered graph
void Node::addHelpNeigh(shared_ptr<Node> n){
    _helpNeigh.push_back(n);
}

// function to remove a neighbour from the list of pointers to the neighbours of a node, NOT TESTED
void Node::removeNeigh(shared_ptr<Node> n){
    _neigh.remove(n);
}

// function to remove all the neighbours of a node
void Node::removeAllNeigh(){
    _neigh.clear();
}

// function to remove all the help-neighbours
void Node::removeAllHelpNeigh(){
    _helpNeigh.clear();
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
        // count the number of opinion 0 and 1 in the opinionlist of the node
        for (int i = 0; i < _opinionlist.size(); i++){
            if (_opinionlist[i] == 0){
                opinion0++;
            }
            else if (_opinionlist[i] == 1){
                opinion1++;
            }
        }


        // count the number of opinions 0 and 1 of the neighbours of the node (only take previously active neighbours into account)
       /* for (shared_ptr<Node> n : _neigh){
            if (n->_wasActive){
                if (n->_opinion == 0){
                    opinion0++; 
                }
                else{ 
                    opinion1++;
                }
            }*/
    }

    if (opinion0 > opinion1){
        if (r >= _resistance){
            _opinion = 0;
        }
        else{
            _opinion = _opinion;
        }
    }
    else if (opinion1 > opinion0){
        if (r >= _resistance){
            _opinion = 1;
        }
        else{
            _opinion = _opinion;
        }
    }
    else{
        _opinion = _opinion;
    }
}

// function that sets the opinion of a node equal to its new opinion
/*void Node::setNewOpinion(){
    _opinion = _newOpinion;
}*/

// function that adds an opinion to the opinionlist
void Node::addOpinion(int opinion){
    _opinionlist.push_back(opinion);
}

// function that sends the opinion of the node to the opinionlist of its neighbours
/*void Node::sendOpinion(){
    for (shared_ptr<Node> n: _neigh){
        n->addOpinion(_opinion);
    }
}*/

// function that deactivates the node
void Node::deactivate(){
    _active = false;
}

// function that sets the activeness of a node
void Node::setActive(bool active){
    _active = active;
}

// function that sets the wasActiveness of a node
/*void Node::setWasActive(){
    _wasActive = _active;
}*/

// function that overwrites the == operator to compare 2 nodes
// ATTENTION: maybe also check opinion and active or not!!
bool operator==(Node n1, Node n2){
    if (n1._index == n2._index){
        return true;
    }
    else{ return false;}
}

// overload << operator to print Node (+ opinion and activeness) to screen
ostream& operator<<(ostream& os, const Node& n){
    os << n._index << '-' << n._opinion << '-' << n._active << ' ';
    return os;
}

