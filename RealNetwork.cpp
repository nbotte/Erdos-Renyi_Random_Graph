// Nina Botte

#include <cmath>
#include <memory>
#include "Graph.h"
#include "Node.h"
#include "Edge.h"
#include "RealNetwork.h"
#include <math.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>
using namespace std;

// implement constructor
Real_World_Network::Real_World_Network(int totalNumberOfNodes, vector<vector<int>> edges){
    _numberOfNodes = totalNumberOfNodes;
    _edges = edges;
    
    _nodelist.resize(_numberOfNodes); 

    makeGraph();
}

void Real_World_Network::makeGraph(){
    // add nodes to nodelist
    for (int i = 0; i < _numberOfNodes; i++){
        // defualt have resistance, opinion and activeness equal to 0;
        double resistance = 0.; // variable that determines the resistance of a node
        int opinion = 0; // variable that determines the opinion of a node
        bool active = 0.; // variable that determines if node is active

        Node n = Node(i, opinion, resistance, active);
        addNode(n);
    }

    // add edges (fill neighlist of each node)
    for (int i = 0; i < _edges.size(); i++){
        int indexInNode = _edges[i][0];
        int indexOutNode = _edges[i][1];
        auto N = make_shared<Node>(_nodelist[indexInNode]);
        auto M = make_shared<Node>(_nodelist[indexOutNode]);
        Edge e = Edge(N, M);
        addEdge(e);
    }
}