// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

#include "Node.h"
#include "Edge.h"
#include <iomanip>
#include <iostream>
#include <list>
#include <vector>
#include <random>
#include <algorithm>
using namespace std;

// attention: try not to copy nodes, becauses no proper copy constructor is written (and I do not see how to do this with neighbours that point to nodes itself)

// default constructor
Node::Node(){};

// implementation of constructor, construct a node by defining its name, opinion, resistance and a list of neighbours (empty at construction)
Node::Node(int index, int opinion, double resistance, bool active){
    _index = index; 
    _neigh = list<int>(); 
    _helpNeigh = list<int>();
    _opinion = opinion; 
    _oldOpinion;
    _resistance = resistance; 
    _active = active;
    _threshold = 0.; // default threshold is 0. (no "stubborness")
    _neighOpinion = vector<int>(); // this is the hidden list (use vector dataset) with the opinions that the neighbors posted since the last time the node was active (see paper 8)
    _cluster = 0; // variable that determines to what cluster the node belongs; mainly used in the SBM, BUT can also be used to optimize modularity etc...; Default zero 

}

// implementation of the getters
int Node::index() const {return _index;}
list<int> Node::neigh() const {return _neigh;}
list<int> Node::helpNeigh() const {return _helpNeigh;}
int Node::opinion() const {return _opinion;}
int Node::oldOpinion() const {return _oldOpinion;}
double Node::resistance() const {return _resistance;}
double Node::threshold() const {return _threshold;}
bool Node::active() const {return _active;}
int Node::cluster() const {return _cluster;}
vector<int> Node::neighOpinion() const {return _neighOpinion;}

// function to add an index of a neighbor to the list of indices of neighbours of a node
void Node::addNeigh(int index){
   _neigh.emplace_back(index);
}

// function to add a help-neighbor when rewiring the clustered graph
void Node::addHelpNeigh(int index){
    _helpNeigh.emplace_back(index);
}

// function to remove an index of a neighbor to the list of indices of neighbors of a node, NOT TESTED
void Node::removeNeigh(int index){
    _neigh.remove(index);
}

// function to remove all the neighbors of a node
void Node::removeAllNeigh(){
    _neigh.clear();
}

// function to remove all the help-neighbors
void Node::removeAllHelpNeigh(){
    _helpNeigh.clear();
}

// function that orders the opinions in the hidden timeline neighOpinion
// here the PR ordering is used (orders according to opinion of node)
void Node::orderOpinionsPR(){
    // first: count the number of opinions equal to own opinion
    int count = 0;
    for (int i = 0; i < _neighOpinion.size(); i++){
        if (_neighOpinion[i] == _opinion){
            count++;
        }
    }

    // then: order the opinions according to your own opinion 
    for (int i = 0; i < count; i++){
        while (_neighOpinion[i] != _opinion){
            std::rotate(_neighOpinion.begin() + i, _neighOpinion.begin() + i + 1, _neighOpinion.end());
        }
    }
}

// function that orders the opinions in the hidden timeline neighOpinion
// here the REC ordering is used (orders according to what is most recent)
void Node::orderOpinionsREC(){
    std::reverse(_neighOpinion.begin(), _neighOpinion.end());
}

// function that orders the opinions in the hidden timeline neighOpinion
// here the REF ordering is used (orders randomly)
void Node::orderOpinionsREF(){
    std::random_shuffle(_neighOpinion.begin(), _neighOpinion.end());
}

// function to change the opinion of the node
void Node::changeOpinion(){
    int opinion0 = 0; // counter that counts the number of neighbours with opinion 0

    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    double o = dis(gen); // generate a random number that will determine the new opinion of the node (see paper 8)
    double r = dis(gen); // generate a random number that will determine if the resistant node changes its opinion or not

    // change the opinion of the active node
    if (_active){
        // first: order the hidden list neighOpinion according to some rule and only look at first 20 opinions
        orderOpinionsREC();

        // count the number of opinion 0 and 1 of the first 20 opinions in the ordered list
        if (_neighOpinion.size() != 0){
            if (_neighOpinion.size() >= 20){
                // only keep max 20 posts
                _neighOpinion.resize(20);
            }
            for (int i = 0; i < _neighOpinion.size(); i++){
                if (_neighOpinion[i] == 0){
                    opinion0++;
                }
            }

            // determine fraction of neighbors with opinion 0 and fraction with opinion 1
            double fraction0 = double(opinion0)/_neighOpinion.size();
            double fraction1 = 1. - fraction0;

            //cout << fraction0 << ' ' << fraction1 << endl;

            // majority model with threshold 
            /*if (fraction0 > fraction1 && fraction0 > _threshold){
                _opinion = 0;
            }
            else if (fraction0 < fraction1 && fraction1 > _threshold){
                _opinion = 1;
            }
            else{
                _opinion = _opinion;
            }*/

            // probabilistic model with all nodes having some fraction of stubbornness
            if (o < fraction0){
                if (r > _resistance){
                    _opinion = 0;
                }
                else{
                    _opinion = _opinion;
                }
            }
            else if (o < fraction0 + fraction1){
                if (r > _resistance){
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
        
    }
}

// function that sets the opinion of a node equal to its new opinion
void Node::setOpinion(int opinion){
    _opinion = opinion;
}

// function that sets the variable oldOpinion 
void Node::setOldOpinion(int opinion){
    _oldOpinion = opinion;
}

// function that sets the resistance of the node
void Node::setResistance(double resistance){
    _resistance = resistance;
}

// function that sets the threshold of the node
void Node::setThreshold(double threshold){
    _threshold = threshold;
}

// function that sets the cluster to which the node belongs
void Node::setCluster(int cluster){
    _cluster = cluster;
}

// function that deactivates the node
void Node::deactivate(){
    _active = false;
}

// function that sets the activeness of a node
void Node::setActive(bool active){
    _active = active;
}

// add opinion to the hidden list of posted opinions of neighbors
void Node::addNeighOpinion(int opinion){
    _neighOpinion.push_back(opinion);
}

// function to remove all the opinions of the neighbors
void Node::removeAllNeighOpinion(){
    _neighOpinion.clear();
}

// function to check if an integer is in a list (used to check if a node has a certain neighbour)
bool Node::containsNeigh(int n){
    return find(_neigh.begin(), _neigh.end(), n) != _neigh.end();
}

// function that overwrites the == operator to compare 2 nodes
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

