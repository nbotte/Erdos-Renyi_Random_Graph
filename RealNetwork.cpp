// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

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
    _numberOfNodes = totalNumberOfNodes; // total number of nodes in the real-world network
    _edges = edges; // edges of the real-world network
    
    // reserve enough memory space
    _nodelist.resize(_numberOfNodes); 

    // construct the real-world network
    makeGraph();
}

void Real_World_Network::makeGraph(){
    // add nodes to nodelist
    for (int i = 0; i < _numberOfNodes; i++){
        // defualt: resistance, opinion and activeness equal to 0;
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

// function that calculates the modularity (returns the calculated modularity)
double Real_World_Network::calculateModularity(vector<vector<int>> communities){
    double mod = 0.; // total modularity
    int L_c = 0; // total number of links in community C
    int k_c = 0; // total degree of nodes in community C
    int L = numberOfEdges(); // total number of links in graph
    // loop over all communities
    for (int i = 0; i < communities.size(); i++){
        // for each community loop over all the constituent nodes
        for (int j = 0; j < communities[i].size(); j++){
            int index = communities[i][j];
            for (int neigh : _nodelist[index].neigh()){
                k_c++; // calculate the total degree of nodes in the community
                // check if neigh is also in the community: if so, add the edge, if not, don't add the edge
                if (find(communities[i].begin(), communities[i].end(), neigh) != communities[i].end()){
                    L_c++;
                }
            }
        }
        mod += (double(L_c)/double(2*L) - pow(double(k_c)/double(2*L), 2));
        L_c = 0;
        k_c = 0;
    }
    return mod;
}
