// Nina Botte

#define _USE_MATH_DEFINES
#include "ErdosRenyi.h"
#include "Graph.h"
#include <cmath>
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <string>
#include <vector>
using namespace std;

#ifndef CLUSTEREDRANDOMER_H
#define CLUSTEREDRANDOMER_H

class Clustered_Random_Network : public Erdos_Renyi_Network{
    int _totalNumberOfNodes; // total number of nodes in the clustered graph
    vector<int> _clusterSizes; // vector that contains the number of nodes for each cluster (sum should equal totalNumberOfNodes)
    vector<double> _edgeProbs; // vector that contains the edge probabilities for each cluster
    double _rewireAddProbability; // probability of rewireing or adding edges between clusters (whether you rewire or add depends on 'type' parameter in constructor)
    string _type; // tells whether you rewire or you add edges to make the clustered graph 

public:
    // define a constructor
    Clustered_Random_Network(int totalNumberOfNodes, vector<int> clusterSizes, vector<double> edgeProbs, double rewireAddProbability, string type);

    // declare member functions of class Clustered_Random_Network 
    vector<int> makeErdosRenyi(int numberOfNodes, double edgeProb, double initOp0Frac, int indexStart, int cluster); // function that makes an ER-graph + returns a vector with the indices of the nodes in that cluster
    void makeGraph(); // function that makes a clustered graph
    void rewireEdges(vector<vector<int>>); // function that rewires the edges of a graph between different clusters, takes a vector of vectors of indices of the nodes of the different clusters as argument
    void addEdges(vector<vector<int>>); // function that adds edges between clusters
    void makeRandomFractionStubborn(double fractionResistant); // function that makes a fraction of randomly selected nodes stubborn
    void makeRandomCommunityFractionStubborn(double fractionResistant); // function that makes a random fraction of communities stubborn
    double calculateModularity(); // function that calculates the modularity of the SBM (based on the clusters that construct the model!) (returns the calculated modularity)

    vector<double> countOpinionFractionCluster(int clusterNumber); // count the fractions of opinions in a particular cluster of the clustered graph
};

#endif