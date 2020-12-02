// Nina Botte

#include <cmath>
#include "ClusteredRandomER.h"
#include "Graph.h"
#include "ErdosRenyi.h"
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <random>
#include <iterator>
#include <string>
#include <functional>
#include <algorithm>
using namespace std;

// seems to be good!

// for now this class seems to work, but can it be implemented nicer/more efficient? what about use of inheritance?
// TO DO: implement opinion dynamics + make random graphs with number of nodes drawn from log-log distribution

Clustered_Random_Network::Clustered_Random_Network(int totalNumberOfNodes, double rewireAddProbability, string type){
    _totalNumberOfNodes = totalNumberOfNodes; // total number of nodes in the clustered graph
    _rewireAddProbability = rewireAddProbability;
    _type = type; // choice between 'rewire' and 'add'

    _nodelist.resize(_totalNumberOfNodes);
    makeGraph();
}

void Clustered_Random_Network::makeGraph(){
    // is there a nicer way to give indexStart a value?
    vector<int> cluster1 = makeErdosRenyi(100, 0.5, 0.5, 0);
    vector<int> cluster2 = makeErdosRenyi(100, 0.5, 0.5, 100);
    vector<int> cluster3 = makeErdosRenyi(100, 0.5, 0.5, 200);
    vector<int> cluster4 = makeErdosRenyi(100, 0.5, 0.5, 300);
    vector<int> cluster5 = makeErdosRenyi(100, 0.5, 0.5, 400);
    vector<int> cluster6 = makeErdosRenyi(100, 0.5, 0.5, 500);
    vector<int> cluster7 = makeErdosRenyi(100, 0.5, 0.5, 600);
    vector<int> cluster8 = makeErdosRenyi(100, 0.5, 0.5, 700);
    vector<int> cluster9 = makeErdosRenyi(100, 0.5, 0.5, 800);
    vector<int> cluster10 = makeErdosRenyi(100, 0.5, 0.5, 900);

    vector<vector<int>> clusters;
    clusters.push_back(cluster1);
    clusters.push_back(cluster2);
    clusters.push_back(cluster3);
    clusters.push_back(cluster4);
    clusters.push_back(cluster5);
    clusters.push_back(cluster6);
    clusters.push_back(cluster7);
    clusters.push_back(cluster8);
    clusters.push_back(cluster9);
    clusters.push_back(cluster10);

    if (_type == "rewire"){
        // rewire edges
        rewireEdges(clusters);
    }
    else if (_type == "add"){
        // add edges between clusters
        addEdges(clusters);
    }
    else{
        cout << "Error: type not equal to one of the two possibilities 'rewire' or 'add'" << endl;
    }
}

// function that makes a clustered graph
// maybe add resistance, opinion, active, etc. as arguments so that you can make clusters with different properties
vector<int> Clustered_Random_Network::makeErdosRenyi(int numberOfNodes, double edgeProb, double initOp0Frac, int indexStart){
    _numberOfNodes = numberOfNodes;
    _edgeProbability = edgeProb;
    _indexStart = indexStart;
    _initOp0Frac = initOp0Frac;
    this->Erdos_Renyi_Network::makeGraph();
    vector<int> cluster;
    for (int i = indexStart; i < indexStart + numberOfNodes; i++){
        cluster.push_back(i);
    }
    return cluster;
}


// function that rewires the edges of the graph with a certain probability --> put this in clustered graph section?
// if you would throw away edgelist, than loop over nodes + loop over each neighbour and for each neighbour change it with some probability (draw random node from nodelist) --> faster run time
void Clustered_Random_Network::rewireEdges(vector<vector<int>> clusters){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    size_t nelem = 1;
    vector<vector<int>> clus; // will contain the random cluster to which you will rewire
    vector<int> out; // will contain the node in that random cluster to which you rewire

    // loop over clusters
    for (int c = 0; c < clusters.size(); c++){
        // move current cluster to beginning
        std::swap(clusters[0], clusters[c]);
        // for each cluster, loop over its nodes
        for (int i = 0; i < clusters[0].size(); i++){
            for (int index : _nodelist[clusters[0][i]].neigh()){
                if (index > _nodelist[clusters[0][i]].index()){
                    double r = dis(gen); // random number that will decide if edge is rewired or not
                    if (r < _rewireAddProbability){
                        // make sure to compile with c++17 (than sample will not give a problem)
                        // choose a random cluster to which you will rewire (current cluster not included)
                        sample(clusters.begin() + 1, clusters.end(), back_inserter(clus), nelem, mt19937{random_device{}()});
                        vector<int> cluster = clus.back(); 
                        // choose a random node whithin that cluster
                        sample(cluster.begin(), cluster.end(), back_inserter(out), nelem, mt19937{random_device{}()});
                        _nodelist[clusters[0][i]].addHelpNeigh(out.back());
                    }
                    else{
                    _nodelist[clusters[0][i]].addHelpNeigh(index);
                    }
                }     
            }
        }
        // move current cluster back to its original position
        std::swap(clusters[0], clusters[c]);
    }

    for (int i = 0; i < _nodelist.size(); i++){
        _nodelist[i].removeAllNeigh();
    }
    for (int i = 0; i < _nodelist.size(); i++){
        for (int index : _nodelist[i].helpNeigh()){
            _nodelist[i].addNeigh(index);
            _nodelist[index].addNeigh(_nodelist[i].index());
        }
        _nodelist[i].removeAllHelpNeigh();
    }
}

void Clustered_Random_Network::addEdges(vector<vector<int>> clusters){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    // loop over clusters
    for (int c = 0; c < clusters.size(); c++){
        // move current cluster to beginning
        std::swap(clusters[0], clusters[c]);
        // for each cluster, loop over its nodes
        for (int i = 0; i < clusters[0].size(); i++){
            // for each node loop over the other clusters, but make sure that you go only once over each possible edge instead of twice (count from c+1 instead of 1)
            for (int k = c+1; k < clusters.size(); k++){
                // loop over the nodes in that cluster
                for (int n = 0; n < clusters[k].size(); n++){
                    double r = dis(gen); // random number that will decide if edge is added between the nodes or not
                    if (r < _rewireAddProbability){
                        // add that node to the neighborlist of the current node and vice versa
                        _nodelist[clusters[0][i]].addNeigh(clusters[k][n]);
                        _nodelist[clusters[k][n]].addNeigh(clusters[0][i]);
                    }
                }
            }
            
        }    
        // move current cluster back to its original position
        std::swap(clusters[0], clusters[c]);
    }
}

// function that checks whether there is an edge between the nodes u and v
bool Clustered_Random_Network::checkEdge(Node u, Node v){
    if (v.containsNeigh(u.index())){
        return true;
    }
    else{ return false;}
}

// function that calculates the local clustering of a node
double Clustered_Random_Network::localClustering(Node n){
    int triangles = 0;
    // check if there is an edge between any pair of neighbours
    for (int index : n.neigh()){
        for (int index2 : n.neigh()){
            if (index2 != index){
                if (checkEdge(_nodelist[index], _nodelist[index2]) || checkEdge(_nodelist[index2], _nodelist[index])){
                    triangles++;
                }
            }
        }
    }
    double localClus = double(triangles)/(double(n.neigh().size())*(double(n.neigh().size())-1));
    return localClus;
}

// function that calculates the average clustering coefficient, brute force for now --> high clustering takes quite a long time to compute, low clustering is fine
double Clustered_Random_Network::averageClustering(){
    double clustering = 0.;
    for (int u = 0; u < _nodelist.size(); u++){
        clustering += localClustering(_nodelist[u]);
    }
    return clustering/_nodelist.size();
}