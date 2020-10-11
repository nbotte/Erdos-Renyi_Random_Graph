// Nina Botte

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

class Node;

class Edge{
    Node* _inNode; // declare inNode variable (pointer to startnode of edge)
    Node* _outNode; // declare outNode variable (pointer to endnode of edge)

public:
    // constructor, construct edge by defining its inNode and its outNode
    Edge(Node* inNode, Node* outNode){_inNode=inNode; _outNode=outNode;}

    // getter, provides access to data member with corresponding name
    Node* const inNode() {return _inNode;}
    Node* const outNode() {return _outNode;}
};

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

// add edgelist again?
class Node{
    int _index; // declare index variable (= name of node)
    list<Node*> _neigh; // declare neigh variable (= vector of pointers to nodes neighbouring the current node)
    //bool _active; // declare active variable (= determines if node is active or not)
    int _opinion; // declare opinion variable (= opinion of node at current time step, choice between 0 and 1)
    int _newOpinion; // declare newOpinion variable (= opinion of node at the next time step, choice between 0 and 1)
    double _resistance; // declare resistance variable (= stubborness of the node, resistance to change his opinion)
    
public:
    // constructor, construct a node by defining its name, opinion, resistance and a vector of neighbours (empty at construction)
    Node(int index, int opinion, double resistance){_index=index; _neigh=list<Node*>(); _opinion=opinion; _newOpinion=opinion; _resistance=resistance;}

    // getter, provides access to data member with corresponding name
    int const index() {return _index;}
    list<Node*> const neigh() {return _neigh;}
    int const opinion() {return _opinion;}
    int const newOpinion() {return _newOpinion;}
    double const resistance() {return _resistance;}

    // function to add a neighbour to the vector of neighbours of a node
    void addNeigh(Node* n){
        _neigh.push_back(n);
    }

    // function to remove a neighbour from the vector of neighbours of a node, NOT TESTED
    void removeNeigh(Node* n){
        _neigh.erase(remove(_neigh.begin(), _neigh.end(), n), _neigh.end());
    }

    // function to change the opinion of the node
    void changeOpinion(){
        int opinion0 = 0; // counter that counts the number of neighbours with opinion 0
        int opinion1 = 0; // counter that counts the number of neighbours with opinion 1

        random_device rd; // will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
        uniform_real_distribution<> dis(0.0, 1.0);

        double r = dis(gen); // generate a random number that will determine if the resistant node will change his opinion or not
 
        // count the number of opinions 0 and 1 of the neighbours of the node
        for (Node* n : _neigh){
            if (n->_opinion == 0){
                opinion0++; 
            }
            else{ 
                opinion1++;
            }
        }
        // change the opinion of the node according to the majority model and if the random number is bigger than the resistance of the node
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

    // function that sets the opinion of a node equal to its new opinion
    void setNewOpinion(){
        _opinion = _newOpinion;
    }
};

// function that overwrites the << operator to print edges to screen
ostream& operator<<(ostream& os, Edge& e){
    os << e.inNode()->index() << '-' << e.outNode()->index() << ' ';
    return os;
}

// function that overwrites the << operator to print nodes (and their opinion) to screen
ostream& operator<<(ostream& os, Node& n){
    os << n.index() << '-' << n.opinion() << ' ';
    return os;
}

class Erdos_Renyi_Network{
    int _numberOfNodes; // total number of nodes in the graph
    double _edgeProbability; // the probablility of having an edge between any pair of nodes
    vector<Node> _nodelist; // vector of nodes in graph
    vector<Edge> _edgelist; // vector of edges in graph

public:
    // constructor, construct graph by making the nodes and the edges with a given probability
    Erdos_Renyi_Network(int numberOfNodes, double edgeProbability){
        _numberOfNodes = numberOfNodes;
        _edgeProbability = edgeProbability;

        // reserve enough memory space for the vectors
        _nodelist.reserve(_numberOfNodes); 
        _edgelist.reserve(pow(_numberOfNodes, 2)); 

        random_device rd; // will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
        uniform_real_distribution<> dis(0.0, 1.0);

        // add nodes to the graph with 50/50 distribution of the 2 possible opinions
        double fractionResistance = 0.5; // set the fraction of stubborn/resistant nodes
        double resistance; // variable that determines the resistance of a node
        int opinion; // variable that determines the opinion of a node
        for (int i = 0; i < _numberOfNodes; i++){
            double k = dis(gen); // random number to determine if node is stubborn
            if (k <= fractionResistance){
                resistance = 1.;
            }
            else{
                resistance = 0.;
            }
            double r = dis(gen); // random number to determine the opinion of a node
            if (r < 0.5){
                opinion = 0;
            }
            else{
                opinion = 1;
            }
            Node n = Node(i, opinion, resistance);
            addNode(n);
        }    
        // add edge between any pair of nodes with a certain probability
        for (int i = 0; i < _nodelist.size(); i++){
            for (int j = i+1; j < _nodelist.size(); j++){
                double r = dis(gen); // draw a random number that will determine whether there is an edge or not
                if (r < _edgeProbability){
                    addEdge(i, j);
                }
            }
        }
    }

    // getter, provides access to data member with corresponding name
    vector<Node> const nodelist() {return _nodelist;}
    vector<Edge> const edgelist() {return _edgelist;}

    // function to add a node to the graph
    void addNode(Node n){
        _nodelist.push_back(n);
    }

    // function to add an edge to the graph
    void addEdge(int indexIn, int indexOut){
        Node* N = &_nodelist[indexIn];
        Node* M = &_nodelist[indexOut];
        Edge e = Edge(N, M);
        // check if edge is already there
        if (contains(_edgelist, e) == false){
            _edgelist.push_back(e);
            _nodelist[indexIn].addNeigh(&_nodelist[indexOut]); // add outNode of edge to neighbours of inNode of edge
            _nodelist[indexOut].addNeigh(&_nodelist[indexIn]); // add inNode of edge to neighbours of outNode of edge

        }
    }

    // function to change the opinions of the nodes in graph based on majority model
    void changeOpinions(){
        for (int i = 0; i < _nodelist.size(); i++){
            _nodelist[i].changeOpinion();
        }
        for (int i = 0; i < _nodelist.size(); i++){
             _nodelist[i].setNewOpinion();
        }
    }

    // function to remove node and edge are not good yet, but are they needed?
    // function to remove an edge from the graph
    // NOT TESTED YET
    /*void removeEdge(Edge* e){
        _edgelist.erase(remove(_edgelist.begin(), _edgelist.end(), e), _edgelist.end()); // remove edge from edgelist
        e->inNode()->removeNeigh(e->outNode()->index()); // remove outNode from the set of neighbours attached to inNode
        e->outNode()->removeNeigh(e->inNode()->index()); // remove inNode from the set of neighbours attached to outNode
    }

    // function to remove a node from the graph
    // NOT TESTED YET
    void removeNode(int index){
        // remove all neighbours from the neighbour list of that node
        _nodelist[index].neigh().clear();
        // remove node from nodelist
        _nodelist.erase(index);
    }*/

    // function that implements a comparison of 2 edges (== operator), returns a bool // PROBABLY NOT NEEDED ANYMORE!
    /*static bool equalEdge(Edge e1, Edge e2){
        if((e1.inNode()->index() == e2.inNode()->index() && e1.outNode()->index() == e2.outNode()->index()) ||
        (e1.inNode()->index() == e2.outNode()->index() && e1.outNode()->index() == e2.inNode()->index())){
            return true;
        }
        else{
            return false;
        }
    }*/

    // function to check if an element is in vector // NEEDS TO BE TESTED!
    bool contains(const vector<Edge> vec, Edge e){
	    return find(vec.begin(), vec.end(), e) != vec.end();
    }

    // function that counts the fraction of the 2 opinions in the graph (returns a vector with the 2 fractions)
    vector<double> countOpinionFraction(){
        int opinion0 = 0;
        int opinion1 = 0;
        vector<double> fractions;
        for (int i = 0; i < _nodelist.size(); i++){
            if (_nodelist[i].opinion() == 0){
                opinion0++;
            }
            else if (_nodelist[i].opinion() == 1){
                opinion1++;
            }
        }
        fractions.push_back(double(opinion0)/double(_numberOfNodes));
        fractions.push_back(double(opinion1)/double(_numberOfNodes));
        return fractions;
    }

    // print function 
    void print(){
        cout << "Nodes: ";
        for (int i = 0; i < _nodelist.size(); ++i){
            cout << _nodelist[i]; 
        }
        cout << endl;
        
        cout << "Edges: ";
        for (int i = 0; i < _edgelist.size(); i++){
            cout << _edgelist[i]; 
        }
        cout << endl;
    }
};

int main(){
    // This example seem to work as expected
   /* Node n = Node(1);
    Node m = Node(2);
    Node k = Node(0);
    Node l = Node(3);

    Edge e1 = Edge(&n, &m);
    Edge e2 = Edge(&n, &k);
    Edge e3 = Edge(&n, &l);
    Edge e4 = Edge(&m, &k);
    Edge e5 = Edge(&l, &k);

    n.addNeigh(m.index());
    m.addNeigh(n.index());
    n.addNeigh(k.index());
    k.addNeigh(n.index());
    n.addNeigh(l.index());
    l.addNeigh(n.index());
    m.addNeigh(k.index());
    k.addNeigh(m.index());
    l.addNeigh(k.index());
    k.addNeigh(l.index());

    n.removeNeigh(m.index());

    vector<int> neigh = n.neigh();
    for (int i = 0; i < neigh.size(); i++){
        cout << neigh[i] << endl;   
    }*/

    // looks ok, but needs further testing
    // possible test: put resistance equal to 1 for every node --> opinion fraction shouldn't change over time
    // N = 1000 and p = 0.1 takes a really long time to run (>1h30min)
    int N = 100;
    double p = 0.1;
    Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p);
    for (int t=0; t<300; t++){
        cout << g.countOpinionFraction()[0] << ' ' << g.countOpinionFraction()[1] << endl;
        g.print();
        cout << endl;
        g.changeOpinions();
    }
    ofstream opfile("Fraction_of_opinions.txt");
    opfile << setprecision(15);
    vector<double> fractionsA(300);
    vector<double> fractionsB(300); 
    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 100; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p); 
        // for each network: let the opinions evolve in time
        for (int t = 0; t < 300; t++){
            double oldFractionA = fractionsA[t];
            double oldFractionB = fractionsB[t];
            fractionsA[t] = oldFractionA + g.countOpinionFraction()[0];
            fractionsB[t] = oldFractionB + g.countOpinionFraction()[1];
            g.changeOpinions();
        }
    }
    // for each time step: print the average opinion fraction over the different graphs to see the opinion evolution
    for (int i = 0; i < 300; i++){
        opfile << fractionsA[i]/100. << ' ' << fractionsB[i]/100. << '\n';
    }
    opfile.close();    
};

// maybe also include adjecency matrix
// https://stackoverflow.com/questions/4964482/how-to-create-two-classes-in-c-which-use-each-other-as-data --> explains 
// how to use 2 classes that reference each other, maybe reconstruct with header files?

// QUESTION: does every class need a destructor? 


/* current status: opinion dynamics: takes often long time to run --> optimizations possible? ALSO: 50/50 case with no stubborn actors almost always leads to a stable 
situation of 54% of one opinion and 46% of the other, however one would expect a 50/50 situation
--> if you take more time steps and average over more simulations, this problem seems to disappear (still needs to be tested in some more depth)*/

// TO DO: play more with fraction of stubborn actors + start implementing active/non-active nodes
// also to do some tests with 1000 nodes
