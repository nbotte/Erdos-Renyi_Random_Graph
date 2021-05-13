// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

#include <cmath>
#include "ClusteredRandom.h"
#include "Graph.h"
#include "ErdosRenyi.h"
#include "WattsStrogatz.h"
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

// constructor for SBM-ER
Clustered_Random_Network::Clustered_Random_Network(int totalNumberOfNodes, vector<int> clusterSizes, vector<double> edgeProbs, double rewireAddProbability, string type){
    _totalNumberOfNodes = totalNumberOfNodes; // total number of nodes in the clustered graph
    _clusterSizes = clusterSizes; // vector that contains the number of nodes for each community (sum should equal totalNumberOfNodes); length of this vector determines number of communities
    _edgeProbs = edgeProbs; // vector that contains the edge probabilities for each community
    _rewireAddProbability = rewireAddProbability; // probability to rewire or add edges between communities; rewiring or adding depends on the 'type' parameter
    _type = type; // choice between 'rewire' and 'add'
    _ER = true; 
    _WS = false;

    // reserve enough memory space
    _nodelist.resize(_totalNumberOfNodes);
    // make a SBM-ER graph
    makeGraph();
}

// constructor for SBM-WS
Clustered_Random_Network::Clustered_Random_Network(int totalNumberOfNodes, vector<int> clusterSizes, vector<double> edgeProbs, vector<int> meanDegrees, double rewireAddProbability, string type){
    _totalNumberOfNodes = totalNumberOfNodes; // total number of nodes in the clustered graph
    _clusterSizes = clusterSizes; // vector that contains the number of nodes for each community (sum should equal totalNumberOfNodes); length of this vector determines number of communities
    _edgeProbs = edgeProbs; // vector that contains the edge probabilities for each community (= the beta parameter for the WS model)
    _meanDegrees = meanDegrees; // vector that contains the mean degrees of each constituent WS model
    _rewireAddProbability = rewireAddProbability; // probability to rewire or add edges between communities; rewiring or adding depends on the 'type' parameter
    _type = type; // choice between 'rewire' and 'add'
    _ER = false;
    _WS = true;

    // reserve enough memory space
    _nodelist.resize(_totalNumberOfNodes);
    // make a SBM-WS graph
    makeGraph();
}

// implementation of getter
vector<vector<int>> Clustered_Random_Network::clusters() const {return _clusters;}

void Clustered_Random_Network::makeGraph(){
    vector<int> cluster;
    int indexStart = 0;
    for (int i = 0; i < _clusterSizes.size(); i++){
        // each community is either an ER graph (for SBM-ER) or a WS graph (for SBM-WS)
        if (_ER == true && _WS == false){
            cluster = makeErdosRenyi(_clusterSizes[i], _edgeProbs[i], 0.5, indexStart, i);
        }
        else if (_ER == false && _WS == true){
            cluster = makeWattsStrogatz(_clusterSizes[i], _meanDegrees[i], _edgeProbs[i], 0.5, indexStart, i);
        }
        _clusters.push_back(cluster); // add the community to the vector that contains the communities
        indexStart += _clusterSizes[i];
    }

    if (_type == "rewire"){
        // rewire edges
        rewireEdges(_clusters);
    }
    else if (_type == "add"){
        // add edges between communities
        addEdges(_clusters);
    }
    else{
        cout << "Error: type not equal to one of the two possibilities 'rewire' or 'add'" << endl;
    }
}

// function that makes the ER communities
vector<int> Clustered_Random_Network::makeErdosRenyi(int numberOfNodes, double edgeProb, double initOp0Frac, int indexStart, int clust){
    _numberOfNodes = numberOfNodes; // number of nodes in specific community
    _edgeProbability = edgeProb; // edge probability of the ER graph that makes the community

    // set indexStart and initial fraction of opinion 0 in the community
    Erdos_Renyi_Network::_indexStart = indexStart;
    Erdos_Renyi_Network::_initOp0Frac = initOp0Frac;

    this->Erdos_Renyi_Network::makeGraph(); // make the community
    vector<int> cluster;
    for (int i = indexStart; i < indexStart + numberOfNodes; i++){
        cluster.push_back(i);
        _nodelist[i].setCluster(clust); // set to which community the nodes belongs
    }
    return cluster; // return the vector that contains the indices of the nodes inside that community
}

// function that makes the WS communities
vector<int> Clustered_Random_Network::makeWattsStrogatz(int numberOfNodes, int meanDegree, double beta, double initOp0Frac, int indexStart, int clust){
    _numberOfNodes = numberOfNodes; // number of nodes in specific community
    _meanDegree = meanDegree; // mean degree in that community
    _rewireProb = beta; // rewire probability of the WS graph that makes the community

    // set indexStart and initial fraction of opinion 0 in the community
    Watts_Strogatz_Network::_indexStart = indexStart;
    Watts_Strogatz_Network::_initOp0Frac = initOp0Frac;

    this->Watts_Strogatz_Network::makeGraph(); // make the community
    vector<int> cluster;
    for (int i = indexStart; i < indexStart + numberOfNodes; i++){
        cluster.push_back(i);
        _nodelist[i].setCluster(clust); // set to which community the node belongs
    }
    return cluster; // return the vector that contains the indices of the nodes inside that community
}

// function that rewires the edges of the graph with a certain probability
void Clustered_Random_Network::rewireEdges(vector<vector<int>> clusters){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    size_t nelem = 1;
    vector<vector<int>> clus; // will contain the random community to which you will rewire
    vector<int> out; // will contain the node in that random community to which you rewire

    // loop over communities
    for (int c = 0; c < clusters.size(); c++){
        // move current community to beginning
        std::swap(clusters[0], clusters[c]);
        // for each community, loop over its nodes
        for (int i = 0; i < clusters[0].size(); i++){
            for (int index : _nodelist[clusters[0][i]].neigh()){
                if (index > _nodelist[clusters[0][i]].index()){
                    double r = dis(gen); // random number that will decide if edge is rewired or not
                    if (r < _rewireAddProbability){
                        // make sure to compile with c++17 (than sample will not give a problem)
                        // choose a random community to which you will rewire (current community not included)
                        sample(clusters.begin() + 1, clusters.end(), back_inserter(clus), nelem, mt19937{random_device{}()});
                        vector<int> cluster = clus.back(); 
                        // choose a random node whithin that community
                        sample(cluster.begin(), cluster.end(), back_inserter(out), nelem, mt19937{random_device{}()});
                        _nodelist[clusters[0][i]].addHelpNeigh(out.back());
                    }
                    else{
                        _nodelist[clusters[0][i]].addHelpNeigh(index);
                    }
                }     
            }
        }
        // move current community back to its original position
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

    // loop over community
    for (int c = 0; c < clusters.size(); c++){
        // move current community to beginning
        std::swap(clusters[0], clusters[c]);
        // for each community, loop over its nodes
        for (int i = 0; i < clusters[0].size(); i++){
            // for each node loop over the other community, but make sure that you go only once over each possible edge instead of twice (count from c+1 instead of 1)
            for (int k = c+1; k < clusters.size(); k++){
                // loop over the nodes in that community
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
        // move current community back to its original position
        std::swap(clusters[0], clusters[c]);
    }
}

// function that makes a random fraction of communities stubborn
void Clustered_Random_Network::makeRandomCommunityFractionStubborn(double fractionResistant){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    int k = 0;
    for (int i = 0; i < _clusterSizes.size(); i++){
        double r = dis(gen); // random number that will determine if community is stubborn or not
        if (r < fractionResistant){
            for (int j = k; j < k + _clusterSizes[i]; j++){
               _nodelist[j].setResistance(1.);
            }
        }
        else{
            for (int j = k; j < k + _clusterSizes[i]; j++){
               _nodelist[j].setResistance(0.);
            }
        }
        k += _clusterSizes[i];
    }
}

// function that gives opinion to nodes in a community (possible to give opinions according to any distribution you want)
// frac0 determines the opinion distribution in that community
// cluster determines which community we are dealing with
// indexStart determines the corresponding position of the nodes in the nodelist
void Clustered_Random_Network::setCommunityOpinion(double frac0, int cluster, int indexStart){
    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    int clusterSize = _clusterSizes[cluster];

    // give nodes in that community their opinion
    for (int i = indexStart; i < indexStart + clusterSize; i++){
        double r = dis(gen); // random number that determines opinion of node
        if (r < frac0){
            _nodelist[i].setOpinion(0);
        }
        else{
            _nodelist[i].setOpinion(1);
        }
    }
}


// function that counts the fractions of opinions in a particular community of the clustered random graph
vector<double> Clustered_Random_Network::countOpinionFractionCluster(int clusterNumber){
    int opinion0 = 0;
    int opinion1 = 0;
    vector<double> fractions;
    int clusterLength = _clusterSizes[clusterNumber];
    // Attention: this for loop only works if all community have equal sizes, if not this should be addapted
    for (int i = (clusterNumber*clusterLength); i < ((clusterNumber + 1)*clusterLength); i++){
        if (_nodelist[i].opinion() == 0){
            opinion0++;
        }
        else if (_nodelist[i].opinion() == 1){
            opinion1++;
        }
    }
    fractions.push_back(double(opinion0)/double(clusterLength));
    fractions.push_back(double(opinion1)/double(clusterLength));
    return fractions;
}

// function that calculates the modularity of the SBM (based on the communities that construct the model!) (returns the calculated modularity)
double Clustered_Random_Network::calculateModularity(){
    int refClus = 0;
    double mod = 0.; // total modularity
    int L_c = 0; // total number of links in community C
    int k_c = 0; // total degree of nodes in community C
    int L = numberOfEdges(); // total number of links in graph
    int i = 0;
    while (i < _nodelist.size()){
        int cluster = _nodelist[i].cluster(); // set current community
        // check if node belongs the the current community
        if (cluster == refClus){
            for (int neigh : _nodelist[i].neigh()){
                k_c++; // calculate the total degree of the nodes in the current community
                // check if neighbor of node belongs to the same community --> if so: add edge to L_c
                if (_nodelist[neigh].cluster() == cluster){
                    L_c++; 
                }
            }
            i++;
        }
        // if node doesn't belong to the current community 
        // --> calculate the modularity of the current community with the obtained values for L_c and k_c
        // --> set L_c and k_c to zero (so you can calculate them for the next community)
        // --> swith to the next community
        else{
            mod += (double(L_c)/double(2*L) - pow(double(k_c)/double(2*L), 2));
            L_c = 0;
            k_c = 0;
            refClus = cluster;
        }
    }
    // this is needed to add the modularity of the last community as well
    mod += (double(L_c)/double(L) - pow(double(k_c)/double(2*L), 2));
    return mod; // return the calculated modularity of the graph
}