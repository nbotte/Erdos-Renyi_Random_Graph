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

// function that performs a community detection, returns the max value of the modularity (modularity of best community division)
// ATTENTION: this works REALLY SLOW --> use python code instead!!
double Real_World_Network::commDetection(){
    vector<vector<int>> communities; // vector that contains the communities (one community is represented as a vector with the indices of the nodes it contains)
    vector<int> comm; // vector that contains one community
    vector<vector<int>> helpComm; // helpvector, will be used to update the community vector

    double mod = 0.; // modularity
    double deltaMod = 0.; // modularity difference after merging 2 communities
    bool start = true; // used to start the while loop

    // first assign each node to a different community
    for (int i = 0; i < _nodelist.size(); i++){
        comm.push_back(i);
        communities.push_back(comm);
        comm.clear();
    }
    mod = calculateModularity(communities);
    int indexCommA; // used to keep track of which communities should be merged
    int indexCommB;
    // as long as the modularity increases, the best community division is not obtained
    while (deltaMod > 0.  || start==true){
        start = false;
        // calculate modularity change for merging each pair of communities, keep largest deltaMod
        for (int i = 0; i < communities.size(); i++){
            for (int k = i+1; k < communities.size(); k++){
                double change = calculateModularityChange(communities[i], communities[k]);
                if (change > deltaMod){
                    deltaMod = change;
                    indexCommA = i;
                    indexCommB = k;
                }
            }
        }
        // merge the communities
        for (int i = 0; i < communities.size(); i++){
            if (i != indexCommA && i != indexCommB){
                helpComm.push_back(communities[i]);
            }
            else if (i == indexCommA){
                for (int j = 0; j < communities[i].size(); j++){
                    comm.push_back(communities[i][j]);
                }
            }
            else if (i == indexCommB){
                for (int j = 0; j < communities[i].size(); j++){
                    comm.push_back(communities[i][j]);
                }
            }
        }
        helpComm.push_back(comm);
        comm.clear();
        communities.clear();
        for (int i = 0; i < helpComm.size(); i++){
            communities.push_back(helpComm[i]);
        }
        helpComm.clear();
        mod += deltaMod;
    }
    cout << calculateModularity(communities) << endl;
    return mod;
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

// funtcion that calculates the modularity change after merging the two communities A and B
double Real_World_Network::calculateModularityChange(vector<int> commA, vector<int> commB){
    int l_AB = 0; // number of edges between A and B
    int k_A = 0; // total degree of nodes in commA
    int k_B = 0; // total degree of nodes in commB
    double deltaMod = 0.; // modularity change after merging A and B
    int L = numberOfEdges(); // total number of edges
    // loop over nodes in community A and check how many edges they have with nodes in community B
    for (int i = 0; i < commA.size(); i++){
        int index = commA[i];
        k_A += _nodelist[index].neigh().size();
        // loop over neighbors and check which neighbor is in commB
        for (int neigh : _nodelist[index].neigh()){
            if (find(commB.begin(), commB.end(), neigh) != commB.end()){
                    l_AB++;
                }
        }
    }
    // calculate k_B
    for (int i = 0; i < commB.size(); i++){
        int index = commB[i];
        k_B += _nodelist[index].neigh().size();
    }
    deltaMod = double(l_AB)/double(L) - (double(k_A)*double(k_B)/(2*pow(double(L), 2)));
    return deltaMod;
}