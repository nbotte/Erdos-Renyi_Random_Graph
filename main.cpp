// Nina Botte

#include <cmath>
//#include "ClusteredRandomER.h"
//#include "ClusteredRandomWS.h"
#include "ClusteredRandom.h"
#include "ErdosRenyi.h"
#include "WattsStrogatz.h"
#include "RealNetwork.h"
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <list>
using namespace std;

// function that calculates the distribution of friends with the same opinion 1 for a particular network (according to paper 8)
void distr_of_friends(){
    int N = 8003;
    int K = 10;
    double p = 0.01;
    double initOp0Frac = 0.5;
    double beta = 0.06; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    // part needed for clustered random network
    /*ifstream file;
    file.open("Community_sizes_powerlaw.txt", ios::in);
    if(!file){
        cout << "No such file";
    }*/
    double p_add = 0.0001;
    /*vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};*/
    vector<int> clusterSizes(53); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(53);
    vector<int> meanDegrees(53);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        if (x < 100){
            edgeProbs.push_back(0.07);
        }
        else{
            edgeProbs.push_back(0.02);
        } 
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 151;
        edgeProbs[i] = 0.025;
        meanDegrees[i] = 4;
    }

    vector<double> fractionAt0(N);
    vector<double> fractionAt500(N);
    // important note: for WS with K = 10, you need to take 11 bins (so that 1 is apart); for other situations 10 bins seems to be better
    double binWidth = 0.1;
    int numberOfBins = 1/binWidth;
    //int numberOfBins = 11;
    vector<int> neighOp1HistAt500(numberOfBins);
    vector<int> neighOp1HistAt0(numberOfBins);

    // needed to calculate variance (see https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/)
    /*double meanAt0op0 = 0.;
    double oldMeanAt0op0 = 0.;
    double meanAt0op1 = 0.;
    double oldMeanAt0op1 = 0.;
    double varAt0op0 = 0.;
    double varAt0op1 = 0.;

    double meanAt500op0 = 0.;
    double oldMeanAt500op0 = 0.;
    double meanAt500op1 = 0.;
    double oldMeanAt500op1 = 0.;
    double varAt500op0 = 0.;
    double varAt500op1 = 0.;*/

    int xAt0op0 = 0;
    int xAt0op1 = 0;
    int xAt500op0 = 0;
    int xAt500op1 = 0;

    vector<int> xVecAt0op0;
    vector<int> xVecAt0op1;
    vector<int> xVecAt500op0;
    vector<int> xVecAt500op1;

    int count = 1;

    // variables that will contain the exact number of nodes with all neighbors having either opinion 0 or opinion 1 (no bin width here)
    int x_0_0 = 0;
    int x_0_1 = 0;
    int x_500_0 = 0;
    int x_500_1 = 0;

    //double mod = 0.;
    //int stubborn = 0;

    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);
    
    ifstream file;
    file.open("lastfm_fin_edgelist.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y){
        edge[0] = x;
        edge[1] = y;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(N, edges);

    // average over different networks
    for (int n = 0; n < 10; n++){
        //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
        //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
        //Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac, 0);
       // Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
        //g.setNodeThreshold(0.);
       // mod += (g.calculateModularity()/10.);
        
        cout << "Graph " << n << endl;

        // for each network average over different simulation
        for (int s = 0; s < 5; s++){
            cout << "Simulation " << s << endl;

            // reset the initial opinions to start a new simulation for the same network + make nodes stubborn
            g.resetInitOpinion(initOp0Frac);
            g.makeRandomFractionStubborn(0.75, 1.); // make all nodes resistant (change how resistant they are from 0, 1)
            //g.makeRandomFractionStubborn(0.1);
            //g.setNodeThreshold(0.);

            // give each community opinions according to predefined distributions
            /*int indexStart = 0;
            for (int i = 0; i < clusterSizes.size(); i++){
                double r = dis(gen); // draw a random number that will determine whether the community has one opinion or not
                if (r < 0.3){
                    g.setCommunityOpinion(1., i, indexStart);
                }
                else{
                    g.setCommunityOpinion(0.286, i, indexStart);
                }
                //g.setCommunityOpinion(0.5, i, indexStart);
                indexStart += clusterSizes[i];
            }*/

           /* for (int i = 0; i < g.nodelist().size(); i++){
                if (g.nodelist()[i].resistance() == 1.){
                    stubborn++;
                }
            }*/

            for (int i = 0; i < g.nodelist().size(); i++){
                int opinion1 = 0;

                if (g.nodelist()[i].neigh().size() != 0){
                    for (int index : g.nodelist()[i].neigh()){
                        if (g.nodelist()[index].opinion() == 1){
                            opinion1++;                    
                        }
                    }
                    fractionAt0[i] = double(opinion1)/double(g.nodelist()[i].neigh().size());
                    // count number of nodes with all neigbors having opinion 0 or opinion 1
                    if (fractionAt0[i] == 0.){
                        x_0_0++;
                        xAt0op0++;
                    }
                    if (fractionAt0[i] == 1.){
                        x_0_1++;
                        xAt0op1++;
                    }
                    
                    // determine position of fraction in the histogram eg 0.111 lies in interval 0-0.1 so position is 0
                    int index0 = int(fractionAt0[i]*10);
                    if (index0 == 10){
                        index0--;
                    }
                    neighOp1HistAt0[index0] += 1;  
                }
            }
            // calculate variance at t = 0 for all neigh having op0 and for all neigh having op1
           /* oldMeanAt0op0 = meanAt0op0;
            meanAt0op0 = meanAt0op0 + (double(xAt0op0) - meanAt0op0)/count;
            varAt0op0 = varAt0op0 + (double(xAt0op0) - meanAt0op0) * (double(xAt0op0) - oldMeanAt0op0);
            oldMeanAt0op1 = meanAt0op1;
            meanAt0op1 = meanAt0op1 + (double(xAt0op1) - meanAt0op1)/count;
            varAt0op1 = varAt0op1 + (double(xAt0op1) - meanAt0op1) * (double(xAt0op1) - oldMeanAt0op1);*/

            xVecAt0op0.push_back(xAt0op0);
            xVecAt0op1.push_back(xAt0op1);

            xAt0op0 = 0;
            xAt0op1 = 0;

            // let opinions evolve in time
            for (int i = 0; i < 500; i++){
                g.setNodesActive(p_bern);
                g.changeOpinions();
                g.deactivateNodes();
            }

            // count fraction of neighbors with opinion 1 at t = 500 + make histogram with number of nodes with a certain fraction of neighbors with opinion 1
            for (int i = 0; i < g.nodelist().size(); i++){
                int opinion1 = 0;

                if (g.nodelist()[i].neigh().size() != 0){
                    for (int index : g.nodelist()[i].neigh()){
                        if (g.nodelist()[index].opinion() == 1){
                            opinion1++;                    
                        }
                    }
                    fractionAt500[i] = double(opinion1)/double(g.nodelist()[i].neigh().size());

                    // count number of nodes with all neigbors having opinion 0 or opinion 1
                    if (fractionAt500[i] == 0.){
                        x_500_0++;
                        xAt500op0++;
                    }
                    if (fractionAt500[i] == 1.){
                        x_500_1++;
                        xAt500op1++;
                    }

                    // determine position of fraction in the histogram eg 0.111 lies in interval 0-0.1 so position is 0
                    int index500 = int(fractionAt500[i]*10);
                    if (index500 == 10){
                        index500--;
                    }
                    neighOp1HistAt500[index500] += 1;
                }
            }
            // calculate variance at t = 500 for all neigh having op0 and for all neigh having op1
            /*oldMeanAt500op0 = meanAt500op0;
            meanAt500op0 = meanAt500op0 + (double(xAt500op0) - meanAt500op0)/count;
            varAt500op0 = varAt500op0 + (double(xAt500op0) - meanAt500op0) * (double(xAt500op0) - oldMeanAt500op0);
            oldMeanAt500op1 = meanAt500op1;
            meanAt500op1 = meanAt500op1 + (double(xAt500op1) - meanAt500op1)/count;
            varAt500op1 = varAt500op1 + (double(xAt500op1) - meanAt500op1) * (double(xAt500op1) - oldMeanAt500op1);*/

            xVecAt500op0.push_back(xAt500op0);
            xVecAt500op1.push_back(xAt500op1);
            
            xAt500op0 = 0;
            xAt500op1 = 0;

            count++; 
        }
    } 
    neighOp1HistAt0[0] = x_0_0;
    neighOp1HistAt0[neighOp1HistAt0.size() - 1] = x_0_1;
    neighOp1HistAt500[0] = x_500_0;
    neighOp1HistAt500[neighOp1HistAt500.size() - 1] = x_500_1;
    ofstream normfile("Hist_500_and_0_fraction_friends_opinion1_real_network_lastfm_REC_8003_fracRes=075_stubb=1_av10x5.txt");
   // ofstream varfile("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_10x100_mean_var_11_bins.txt");
    ofstream xfile("Hist_500_and_0_fraction_friends_opinion1_real_network_lastfm_REC_8003_fracRes=075_stubb=1_av10x5_xvalues.txt");
    ofstream echofile("Echo_chamber_real_network_lastfm_REC_8003_fracRes=075_stubb=1_av10x5.txt");
    for (int i = 0; i < neighOp1HistAt500.size(); i++){
        double norm = double(neighOp1HistAt500[i]) / double(neighOp1HistAt0[i]);
        normfile << neighOp1HistAt500[i] << ' ' << neighOp1HistAt0[i] << ' ' << norm << endl;
    }
    //varfile << meanAt0op0 << ' ' << varAt0op0/double(count - 2) << '\n' << meanAt0op1 << ' ' << varAt0op1/double(count - 2) << '\n' <<  meanAt500op0 << ' ' << varAt500op0/double(count - 2) << '\n' << meanAt500op1 << ' ' << varAt500op1/double(count - 2) << endl;
    for (int i = 0; i < 100; i++){
        xfile << xVecAt0op0[i] << ' ' << xVecAt0op1[i] << ' ' << xVecAt500op0[i] << ' ' << xVecAt500op1[i] << endl;
    }
    double echo_0 = double(x_500_0)/double(x_0_0);
    double echo_1 = double(x_500_1)/double(x_0_1);
    echofile << x_0_0 << ' ' << x_500_0 << ' ' << echo_0 << '\n' << x_0_1 << ' ' << x_500_1 << ' ' << echo_1 << endl;
    normfile.close();
    //varfile.close();
    xfile.close();
    echofile.close();
    /*stubborn = stubborn/100;
    cout << stubborn << endl;*/

    //cout << mod << endl;
}

// function that calculates the evolution of opinions in a particular network
void evolution_of_opinions(){
    int N = 8003;
    int K = 10;
    double p = 0.01;
    double initOp0Frac = 0.5;
    double beta = 0.06; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    // part needed for clustered random network
    /*ifstream file;
    file.open("Community_sizes_powerlaw.txt", ios::in);
    if(!file){
        cout << "No such file";
    }*/
    double p_add = 0.0001;
    /*vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};*/
    vector<int> clusterSizes(53); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(53);
    vector<int> meanDegrees(53);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        if (x < 100){
            edgeProbs.push_back(0.07);
        }
        else{
            edgeProbs.push_back(0.02);
        }
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 151;
        edgeProbs[i] = 0.025;
        meanDegrees[i] = 4;
    }
  
    ofstream opfile("Fraction_of_opinions_active_01_av_good_init_real_network_lastfm_REC_8003_fracRes=075_stubb=1_av10x5.txt");
    vector<double> mean0(500); // contains the average fraction of opinion 0 in the graph at each timestep
    vector<double> mean1(500); // contains the average fraction of opinion 1 in the graph at each timestep
    vector<double> variance0(500); // calculate variance of opinion 0 according to Welford's algorithm
    vector<double> variance1(500); // calculate variance of opinion 1 according to Welford's algorithm
    int count = 1;

   // vector<double> commOp0Begin(clusterSizes.size()); // vector that contains the fraction of opinion 0 in each community at beginning of simulation --> also check for n, s = 1 (to see things without averaging)
    //vector<double> commOp0End(clusterSizes.size()); // vector that contains the fraction of opinion 0 in each community at end of simulation --> also check for n, s = 1 (to see things without averaging)

    random_device rd; // will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // standard mersenne twister engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);

    ifstream file;
    file.open("lastfm_fin_edgelist.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y){
        edge[0] = x;
        edge[1] = y;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(N, edges);

    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 1; n++){
       // Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
       // Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
        //Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac, 0);
        //Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
       // g.setNodeThreshold(0.);
        cout << "Graph: " << n << endl;
        // for each network, run different simulations --> can this be implemented faster?
        for (int s = 0; s < 1; s++){
            // reset the initial opinions
            g.resetInitOpinion(initOp0Frac);
          //  g.setNodeThreshold(0.);
            g.makeRandomFractionStubborn(0.75, 1.); // make all nodes resistant (change how resistant they are from 0, 1)

            // give each community opinions according to predefined distributions
            /*int indexStart = 0;
            for (int i = 0; i < clusterSizes.size(); i++){
                double r = dis(gen); // draw a random number that will determine whether the community has one opinion or not
                if (r < 0.3){
                    g.setCommunityOpinion(1., i, indexStart);
                }
                else{
                    g.setCommunityOpinion(0.286, i, indexStart);
                }
                //g.setCommunityOpinion(0.5, i, indexStart);
                indexStart += clusterSizes[i];
            }*/

            // At beginning of each time evolution: calculate fraction of opinion 0 in each cluster
          /*  for (int i = 0; i < clusterSizes.size(); i++){
                double frac0 = g.countOpinionFractionCluster(i)[0];
                commOp0Begin[i] += frac0;
            }*/

            // for each network and each simulation: let the opinions evolve in time
            for (int t = 0; t < 500; t++){
                g.setNodesActive(p_bern);
                g.changeOpinions();

                double x0 = g.countOpinionFraction()[0];
                double x1 = g.countOpinionFraction()[1];

                double oldMean0 = mean0[t];
                double oldMean1 = mean1[t];
                mean0[t] = mean0[t] + (x0 - mean0[t])/count;
                mean1[t] = mean1[t] + (x1 - mean1[t])/count;

                variance0[t] = variance0[t] + (x0 - mean0[t]) * (x0 - oldMean0);
                variance1[t] = variance1[t] + (x1 - mean1[t]) * (x1 - oldMean1);

                g.deactivateNodes();
            }
            // At end of each time evolution: calculate fraction of opinion 0 in each cluster
           /* for (int i = 0; i < clusterSizes.size(); i++){
                double frac0 = g.countOpinionFractionCluster(i)[0];
                commOp0End[i] += frac0;
            }*/
            count++;
        }
    }
    // for each time step: print the average opinion fraction over the different graphs to see the opinion evolution
    for (int i = 0; i < 500; i++){
        opfile << mean0[i] << ' ' << mean1[i] << ' ' << variance0[i]/double(count - 2) << ' ' << variance1[i]/double(count - 2) << endl;
    }
    opfile.close();

    // write community opinion 0 fractions at end of time evolution to file
   /* ofstream commfile("Fraction_of_opinion0_comm_SBM-WS_active_01_av_good_init_REC_10-001-0001_10x100_T=0.txt");
    for (int i = 0; i < commOp0End.size(); i++){
        double frac0Begin = commOp0Begin[i]/100.;
        double frac0End = commOp0End[i]/100.;
        commfile << frac0Begin << ' ' << frac0End << endl;
    }
    commfile.close();*/
}

// function that calculates the degree distribution + the average degree distribution of a particular network
void degree_distr(){
    int N = 8003;
    int K = 6;
    double p = 0.01;
    double initOp0Frac = 0.5;
    double beta = 0.01; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    // part needed for clustered random network
   /* double p_add = 0.001;
    vector<int> clusterSizes(10); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(10);
    for (int i = 0; i < 10; i++){
        clusterSizes[i] = 100;
        edgeProbs[i] = 0.1;
    }*/
    
    ifstream file;
    file.open("lastfm_fin_edgelist.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y){
        edge[0] = x;
        edge[1] = y;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(N, edges);

    vector<int> degreeDistr(70);
   // Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
    for (int i = 0; i < g.nodelist().size(); i++){
        degreeDistr[g.nodelist()[i].neigh().size()] += 1;
    }
    ofstream degreeFile("Degree_distribution_real_network_lastfm_fin.txt");
    for (int i = 0; i < degreeDistr.size(); i++){
        degreeFile << i << ' ' << degreeDistr[i] << '\n';
    }

    degreeFile.close();
  /*  vector<int> degreeDistrAv(10000);
    for (int k = 0; k < 100; k++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
        cout << "Graph: " << k << endl;
        for (int i = 0; i < g.nodelist().size(); i++){
            degreeDistrAv[g.nodelist()[i].neigh().size()] += 1;
        }
    }


    ofstream degreeFileAv("Degree_distribution_real_network_facebook_av.txt");
    for (int i = 0; i < degreeDistrAv.size(); i++){
        degreeFileAv << i << ' ' << degreeDistrAv[i]/100 << '\n';
    }

    degreeFileAv.close();*/
}

// function that calculates the average degree of a graph
void Av_degree(){
    int N = 1000;
    int K = 6;
    double p = 0.01;
    double initOp0Frac = 0.5;
    double beta = 0.01; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    // part needed for clustered random network
    /*ifstream file;
    file.open("Community_sizes_powerlaw.txt", ios::in);
    if(!file){
        cout << "No such file";
    }*/
    double p_add = 0.008;
   /* vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};*/
    vector<int> clusterSizes(10); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(10);
    vector<int> meanDegrees(10);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        if (x < 100){
            edgeProbs.push_back(0.07);
        }
        else{
            edgeProbs.push_back(0.02);
        }        
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 100;
        edgeProbs[i] = 0.03;
        meanDegrees[i] = 10;
    }
   /* ifstream file;
    file.open("lastfm_fin_edgelist.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y){
        edge[0] = x;
        edge[1] = y;
        edges.push_back(edge);
    }
    file.close();*/

    //Real_World_Network g = Real_World_Network(8003, edges);

    //int degree = 0;
    //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
   // Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
    /*for (int i = 0; i < g.nodelist().size(); i++){
        degree += g.nodelist()[i].neigh().size();
    }
    double mod = g.calculateModularity();
    double modTest = g.calculateModularityTest(g.clusters());
    double Av_degree = double(degree)/double(N);*/
    ofstream clusFile("avDeg_SBM_low.txt");
    for (int i = 0; i < 100; i++){
        //Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac, 0);
        Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
        //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
        //Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
        //double mod = g.calculateModularity();
        int degree = 0;
        for (int i = 0; i < g.nodelist().size(); i++){
            degree += g.nodelist()[i].neigh().size();
        }
        double Av_degree = double(degree)/double(N);
        clusFile << Av_degree << '\n';

    }
    clusFile.close();
    
    
   /* int edges_tot = g.numberOfEdges();
    cout << clus << ' ' << Av_degree << ' ' << mod << ' ' << modTest << ' ' << edges_tot << endl;*/
}

// function to do small tests
void test(){
    int N = 100;
    int K = 4;
    double beta = 0.05;
    double initOp0Frac = 0.5;
    double p_bern = 1.;

    // part needed for clustered random network
    /*ifstream file;
    file.open("Community_sizes_powerlaw.txt", ios::in);
    if(!file){
        cout << "No such file";
    }*/
    double p_add = 0.001;
    vector<int> clusterSizes(10); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(10);
    vector<int> meanDegrees(10);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        edgeProbs.push_back(0.5);
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 10;
        edgeProbs[i] = 0.01;
        meanDegrees[i] = 4;
    }

    ifstream file;
    file.open("facebook_data.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y){
        edge[0] = x;
        edge[1] = y;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(4039, edges);
    //cout << g.commDetection() << endl;

    //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
    //Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac, 0);
   //Erdos_Renyi_Network g = Erdos_Renyi_Network(N, beta, p_bern, initOp0Frac, 0);
   // cout << g.calculateModularity() << ' ' << g.averageClustering() << endl;
   /* int k = 0;
    for (int i = 0; i < clusterSizes.size(); i++){
        if (i%2){
            g.setCommunityOpinion(0.5, i, k);
        }
        else{
           g.setCommunityOpinion(0.5, i, k); 
        }
        k += clusterSizes[i];
    }
    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ' ';
    }
    cout << endl;*/


    /*g.makeRandomFractionStubborn(0.3);
    int stubborn = 0;
    for (int i = 0; i < g.nodelist().size(); i++){
        if (g.nodelist()[i].resistance() == 1.){
            stubborn++;
        }
    }
    cout << stubborn << endl;*/
    //cout << g.calculateModularity() << endl;
   // g.setNodesActive(p_bern);
    //g.setNodeThreshold(0.3);

   /* for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ' ' << g.nodelist()[i].cluster() << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index] << '\t';
        }
        cout << endl;
    }*/

   /* g.changeOpinions();
    g.changeOpinions();

    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index] << '\t';
        }
        cout << endl;
    }*/

    vector<vector<int>> comm = {{3304, 1272, 3328, 3330, 3331, 3336, 3338, 3341, 3349, 3353, 3358, 1325, 1327, 3366, 3372, 1339, 3375, 3377, 2661, 2662, 2665, 2667, 2668, 2671, 2672, 2674, 2675, 2677, 2678, 2682, 2684, 2685, 2686, 2687, 2688, 3383, 2690, 2691, 2692, 2693, 2696, 2697, 2700, 2701, 2702, 2707, 2708, 2709, 2711, 2712, 2713, 2714, 2715, 2716, 2718, 2719, 2720, 2721, 2722, 2723, 2724, 2725, 2727, 2728, 1017, 2730, 2732, 2733, 2734, 2735, 2736, 2737, 2739, 2742, 2744, 2747, 2748, 2751, 2752, 2753, 3396, 3397, 2758, 2759, 2760, 2762, 2764, 2765, 2768, 2770, 2771, 2772, 2775, 2776, 2778, 2779, 2780, 2781, 2783, 2784, 2785, 2787, 2788, 2790, 2791, 2792, 2793, 2795, 2796, 2797, 2799, 2801, 2802, 2803, 2804, 2805, 2808, 2811, 2812, 2813, 2816, 2818, 2819, 2820, 2822, 2823, 2824, 2825, 2826, 2828, 2829, 2830, 2831, 2832, 2833, 2836, 2840, 2841, 2842, 2843, 2844, 2845, 2846, 2847, 2848, 2849, 3415, 2853, 2855, 2856, 2857, 2858, 2859, 2860, 2861, 3418, 2863, 2865, 2868, 2869, 2870, 2871, 2875, 2876, 2878, 2881, 2882, 2883, 2886, 2889, 2892, 2893, 2894, 2895, 2898, 2900, 2901, 2902, 2903, 2904, 2914, 2918, 2921, 2922, 2923, 2926, 2930, 2932, 2933, 2934, 2935, 2936, 2937, 2938, 2939, 2941, 2942, 897, 2946, 899, 2947, 2949, 2950, 2951, 2952, 906, 2954, 2955, 908, 2957, 2958, 2961, 2962, 916, 2964, 2965, 920, 921, 2968, 922, 2970, 925, 926, 923, 2971, 2973, 927, 2975, 928, 932, 2976, 929, 934, 2979, 2980, 2983, 2984, 2989, 2990, 2991, 2992, 946, 947, 2995, 2996, 2997, 950, 952, 953, 2998, 951, 2999, 3001, 3002, 959, 3005, 3006, 960, 3007, 961, 3009, 3011, 3012, 3013, 966, 967, 3015, 972, 3016, 3017, 970, 3018, 3019, 3020, 973, 3021, 3024, 978, 3031, 3027, 982, 983, 980, 3028, 3032, 3034, 3035, 3037, 990, 991, 3038, 993, 3039, 995, 998, 3041, 3043, 997, 996, 3046, 999, 3048, 3049, 3050, 1003, 3051, 1004, 1006, 3053, 1008, 3063, 3057, 3058, 3066, 3067, 3060, 3069, 3061, 1023, 1024, 3072, 1026, 3075, 1028, 1029, 3076, 3077, 1031, 3079, 1033, 1034, 3081, 3082, 3083, 1039, 1038, 1040, 3086, 3087, 3088, 3089, 3091, 3093, 3094, 1047, 1049, 1048, 3100, 3101, 3102, 3103, 1054, 3105, 3106, 1059, 3107, 1061, 3109, 3110, 1063, 3113, 3111, 3112, 1068, 1069, 3062, 3118, 3115, 3121, 3122, 3116, 3117, 1074, 1077, 3120, 1073, 1075, 1076, 1078, 1084, 1079, 1083, 3128, 3131, 3132, 1086, 1091, 1092, 3134, 1087, 1088, 3138, 3140, 3141, 3144, 3068, 3145, 1098, 3146, 3148, 1101, 3150, 3151, 3152, 3154, 3070, 1110, 1112, 1107, 3155, 3071, 3156, 1116, 3157, 3159, 3160, 3161, 3162, 3164, 1117, 3165, 3167, 3168, 3169, 3170, 1123, 1124, 1125, 1133, 1126, 3183, 1128, 3176, 1132, 1135, 3179, 3180, 3181, 3182, 1136, 3185, 3186, 3187, 1140, 3188, 3189, 3190, 1144, 1146, 1153, 1149, 1156, 3196, 3197, 1150, 1160, 3200, 3201, 1163, 1164, 3204, 1157, 3206, 3207, 3208, 1161, 3210, 1172, 1173, 1165, 3214, 3215, 3217, 3218, 3221, 1180, 1181, 3230, 1175, 3224, 1184, 3227, 1050, 1185, 1182, 3233, 1191, 1187, 3235, 3236, 3239, 3240, 3241, 3242, 1198, 1195, 1201, 1196, 3244, 1199, 3245, 1205, 1207, 3251, 1209, 3258, 3259, 3260, 1055, 3261, 3262, 3263, 3264, 1056, 3265, 3266, 3267, 3270, 3268, 1220, 3271, 1222, 3273, 3272, 3276, 3274, 3279, 1058, 3278, 3275, 3277, 3282, 1230, 1238, 1239, 3281, 3280, 3285, 3287, 3289, 3288, 1243, 3097, 3296, 3297, 3291, 1242, 1250, 3292, 3298, 3295, 1251, 3299, 3306, 3305, 1256, 3307, 3302, 3099, 1255, 3313, 3308, 3312, 1266, 1265, 3310, 1271, 1267, 3320, 3316, 1274, 1269, 3319, 3326, 3322, 3323, 3321, 3329, 3327, 1283, 1285, 1280, 1278, 3334, 3337, 1290, 3335, 3340, 3332, 1287, 1288, 3342, 3343, 1289, 3347, 1291, 1301, 1293, 3344, 3345, 3346, 3350, 3355, 3354, 3351, 3352, 1302, 3360, 3359, 1305, 1307, 3356, 3365, 3357, 3361, 3368, 1312, 3370, 1323, 1317, 3367, 1321, 3369, 3371, 1329, 1330, 3373, 3374, 1328, 3376, 3378, 3379, 1331, 3386, 1335, 1336, 3380, 3381, 3385, 1341, 1337, 1346, 1344, 1340, 3389, 3391, 1351, 3393, 3395, 3402, 3403, 3398, 3405, 3401, 1352, 1359, 3408, 3409, 3407, 3404, 1361, 3410, 3406, 3412, 1360, 1367, 1370, 3411, 1365, 3422, 1369, 3416, 3420, 3421, 3425, 3419, 1375, 3427, 3424, 1376, 1377, 3124, 1380, 3434, 1383, 3429, 3432, 3431, 3125, 1388, 1389, 3126, 1393, 1390, 1391, 1398, 1401, 1399, 1402, 1405, 1407, 1409, 1410, 1411, 1416, 1419, 1420, 1421, 3133, 1431, 1433, 1434, 1437, 1439, 1440, 1441, 1442, 3136, 1445, 1447, 1449, 1450, 1456, 1457, 1458, 1460, 1461, 1463, 1467, 1470, 1471, 1476, 1477, 1480, 1483, 1484, 1485, 1488, 1491, 1494, 1498, 1501, 1505, 1509, 1511, 1513, 1516, 1517, 1518, 1519, 1520, 1521, 1522, 1523, 1524, 1527, 1528, 1530, 1532, 1534, 1535, 1537, 1538, 1539, 1542, 1544, 1547, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1559, 1561, 1563, 1564, 1567, 1570, 1571, 1572, 1579, 1580, 1584, 1585, 1588, 1589, 1130, 1590, 1593, 1594, 1597, 1598, 1600, 1603, 1604, 1605, 1608, 1609, 1610, 1612, 1613, 1614, 1615, 1617, 1618, 3171, 1619, 1620, 1621, 1622, 1623, 3172, 1626, 3173, 1632, 1637, 3175, 1639, 1641, 1642, 1643, 1644, 3177, 1651, 1652, 1653, 3178, 1656, 1659, 1662, 1663, 1665, 1666, 1668, 1669, 1670, 1674, 1675, 1676, 1678, 1683, 1684, 1685, 1687, 1688, 1689, 1697, 1698, 1700, 1701, 1702, 1705, 1707, 1708, 1710, 1712, 1714, 1717, 3191, 1719, 1721, 1722, 1723, 1724, 3192, 1726, 1730, 1734, 3194, 1735, 1736, 1737, 1741, 1746, 1750, 1752, 1753, 1754, 3198, 1757, 1758, 3199, 1761, 1765, 1768, 1769, 1771, 1772, 1774, 1775, 1779, 3203, 1780, 1782, 1789, 3205, 1791, 1792, 1793, 1795, 1796, 1797, 1799, 1800, 1803, 1804, 1805, 1806, 1809, 3209, 1811, 1810, 1174, 1813, 1816, 1817, 1819, 3211, 1821, 1822, 1823, 3212, 1825, 1826, 1827, 1178, 1832, 1833, 1835, 1836, 1838, 1839, 1842, 1843, 3216, 1845, 1846, 1849, 1851, 1852, 1854, 1856, 1858, 3219, 1860, 1861, 1863, 1864, 1865, 1866, 1867, 1868, 1874, 3222, 1877, 1879, 3223, 1883, 1886, 1888, 3225, 1891, 1898, 1900, 1902, 1905, 1908, 1909, 1911, 3246, 1211, 3247, 3248, 3249, 1214, 3252, 3253, 3254, 1219, 3255, 3256}, 
    {107, 348, 349, 350, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 394, 395, 396, 397, 398, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 455, 456, 457, 458, 459, 460, 461, 462, 463, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 579, 580, 581, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 596, 597, 598, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 629, 630, 631, 633, 634, 636, 637, 638, 639, 641, 642, 644, 645, 646, 648, 649, 651, 652, 653, 654, 655, 656, 657, 660, 663, 664, 666, 667, 668, 669, 671, 672, 673, 674, 676, 677, 678, 679, 680, 682, 683, 684, 685, 896, 898, 900, 902, 904, 905, 907, 909, 910, 911, 912, 913, 914, 915, 917, 918, 919, 924, 930, 931, 933, 935, 936, 937, 939, 940, 941, 942, 943, 944, 945, 948, 949, 954, 955, 956, 957, 962, 964, 965, 968, 969, 971, 974, 975, 976, 977, 979, 981, 984, 985, 986, 987, 988, 989, 994, 1000, 1001, 1002, 1005, 1007, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1018, 1019, 1020, 1021, 1022, 1025, 1027, 1030, 1032, 1035, 1036, 1037, 1041, 1042, 1043, 1044, 1045, 1046, 1051, 1052, 1053, 1057, 1060, 1062, 1064, 1065, 1066, 1067, 1070, 1071, 1072, 1080, 1081, 1082, 1089, 1090, 1093, 1094, 1095, 1096, 1099, 1100, 1102, 1103, 1104, 1105, 1106, 1108, 1109, 1111, 1113, 1114, 1115, 1119, 1120, 1121, 1122, 1127, 1129, 1131, 1134, 1138, 1139, 1141, 1142, 1143, 1145, 1147, 1148, 1151, 1152, 1154, 1155, 1158, 1159, 1162, 1166, 1167, 1168, 1169, 1170, 1176, 1177, 1179, 1183, 1188, 1189, 1190, 1192, 1194, 1197, 1200, 1202, 1203, 1204, 1206, 1208, 1210, 1212, 1213, 1215, 1217, 1218, 1221, 1223, 1224, 1225, 1226, 1227, 1228, 1229, 1231, 1233, 1234, 1235, 1236, 1237, 1240, 1241, 1244, 1245, 1246, 1247, 1248, 1249, 1252, 1253, 1254, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264, 1268, 1270, 1273, 1275, 1276, 1277, 1279, 1281, 1282, 1284, 1286, 1292, 1294, 1295, 1296, 1298, 1299, 1300, 1303, 1304, 1306, 1308, 1309, 1310, 1311, 1313, 1315, 1316, 1318, 1320, 1322, 1324, 1326, 1332, 1334, 1338, 1342, 1343, 1345, 1347, 1348, 1349, 1350, 1353, 1354, 1355, 1356, 1357, 1358, 1362, 1364, 1366, 1368, 1371, 1372, 1373, 1374, 1379, 1381, 1382, 1384, 1385, 1386, 1392, 1394, 1395, 1396, 1397, 1400, 1403, 1404, 1406, 1408, 1412, 1413, 1414, 1415, 1417, 1418, 1422, 1423, 1425, 1426, 1427, 1428, 1429, 1430, 1432, 1435, 1436, 1438, 1443, 1444, 1446, 1448, 1451, 1453, 1454, 1455, 1459, 1462, 1464, 1466, 1469, 1472, 1473, 1474, 1475, 1478, 1479, 1481, 1482, 1487, 1489, 1490, 1492, 1495, 1496, 1497, 1499, 1500, 1502, 1503, 1504, 1506, 1507, 1508, 1510, 1512, 1514, 1515, 1525, 1526, 1529, 1531, 1536, 1540, 1541, 1543, 1545, 1546, 1550, 1558, 1560, 1562, 1565, 1566, 1569, 1573, 1574, 1575, 1576, 1578, 1581, 1582, 1583, 1586, 1587, 1591, 1592, 1595, 1596, 1599, 1602, 1606, 1607, 1611, 1616, 1624, 1625, 1627, 1628, 1630, 1631, 1633, 1634, 1635, 1636, 1638, 1640, 1645, 1646, 1647, 1649, 1650, 1654, 1655, 1657, 1658, 1660, 1661, 1664, 1667, 1671, 1672, 1673, 1677, 1679, 1680, 1681, 1682, 1686, 1690, 1691, 1692, 1693, 1694, 1695, 1696, 1699, 1703, 1704, 1706, 1709, 1711, 1713, 1715, 1716, 1720, 1725, 1727, 1728, 1729, 1731, 1732, 1738, 1739, 1740, 1742, 1743, 1744, 1745, 1748, 1749, 1751, 1755, 1756, 1759, 1760, 1762, 1763, 1764, 1766, 1767, 1770, 1773, 1776, 1777, 1778, 1781, 1783, 1785, 1786, 1788, 1790, 1794, 1801, 1807, 1808, 1812, 1814, 1815, 1818, 1820, 1824, 1828, 1829, 1830, 1831, 1834, 1840, 1841, 1844, 1847, 1848, 1850, 1853, 1855, 1857, 1859, 1862, 1869, 1870, 1871, 1872, 1873, 1875, 1876, 1878, 1881, 1882, 1884, 1885, 1887, 1889, 1890, 1893, 1894, 1896, 1897, 1899, 1901, 1903, 1904, 1906, 1907, 1910, 1967},
    {857, 862, 865, 868, 1085, 3437, 3438, 3439, 3440, 3441, 3442, 3443, 3444, 3445, 3446, 3447, 3448, 3449, 3450, 3451, 3452, 3453, 3454, 3455, 3456, 3457, 3458, 3459, 3460, 3461, 3462, 3463, 3464, 3465, 3466, 3467, 3468, 3469, 3470, 3471, 3472, 3473, 3474, 3475, 3476, 3477, 3478, 3479, 3480, 3481, 3482, 3483, 3484, 3485, 3486, 3487, 3488, 3489, 3490, 3491, 3492, 3493, 3494, 3495, 3496, 3497, 3498, 3499, 3500, 3501, 3502, 3503, 3504, 3505, 3506, 3507, 3508, 3509, 3510, 3511, 3512, 3513, 3514, 3515, 3516, 3517, 3518, 3519, 3520, 3521, 3522, 3523, 3524, 3525, 3526, 3527, 3528, 3529, 3530, 3531, 3532, 3533, 3534, 3535, 3536, 3537, 3538, 3539, 3540, 3541, 3542, 3543, 3544, 3545, 3546, 3547, 3548, 3549, 3550, 3551, 3552, 3553, 3554, 3555, 3556, 3557, 3558, 3559, 3560, 3561, 3562, 3563, 3564, 3565, 3566, 3567, 3568, 3569, 3570, 3571, 3572, 3573, 3574, 3575, 3576, 3577, 3578, 3579, 3580, 3581, 3582, 3583, 3584, 3585, 3586, 3587, 3588, 3589, 3590, 3591, 3592, 3593, 3594, 3595, 3596, 3597, 3598, 3599, 3600, 3601, 3602, 3603, 3604, 3605, 3606, 3607, 3608, 3609, 3610, 3611, 3612, 3613, 3614, 3615, 3616, 3617, 3618, 3619, 3620, 3621, 3622, 3623, 3624, 3625, 3626, 3627, 3628, 3629, 3630, 3631, 3632, 3633, 3634, 3635, 3636, 3637, 3638, 3639, 3640, 3641, 3642, 3643, 3644, 3645, 3646, 3647, 3648, 3649, 3650, 3651, 3652, 3653, 3654, 3655, 3656, 3657, 3658, 3659, 3660, 3661, 3662, 3663, 3664, 3665, 3666, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3674, 3675, 3676, 3677, 3678, 3679, 3680, 3681, 3682, 3683, 3684, 3685, 3686, 3687, 3688, 3689, 3690, 3691, 3692, 3693, 3694, 3695, 3696, 3697, 3698, 3699, 3700, 3701, 3702, 3703, 3704, 3705, 3706, 3707, 3708, 3709, 3710, 3711, 3712, 3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725, 3726, 3727, 3728, 3729, 3730, 3731, 3732, 3733, 3734, 3735, 3736, 3737, 3738, 3739, 3740, 3741, 3742, 3743, 3744, 3745, 3746, 3747, 3748, 3749, 3750, 3751, 3752, 3753, 3754, 3755, 3756, 3757, 3758, 3759, 3760, 3761, 3762, 3763, 3764, 3765, 3766, 3767, 3768, 3769, 3770, 3771, 3772, 3773, 3774, 3775, 3776, 3777, 3778, 3779, 3780, 3781, 3782, 3783, 3784, 3785, 3786, 3787, 3788, 3789, 3790, 3791, 3792, 3793, 3794, 3795, 3796, 3797, 3798, 3799, 3800, 3801, 3802, 3803, 3804, 3805, 3806, 3807, 3808, 3809, 3810, 3811, 3812, 3813, 3814, 3815, 3816, 3817, 3818, 3819, 3820, 3821, 3822, 3823, 3824, 3825, 3826, 3827, 3828, 3829, 3830, 3831, 3832, 3833, 3834, 3835, 3836, 3837, 3838, 3839, 3840, 3841, 3842, 3843, 3844, 3845, 3846, 3847, 3848, 3849, 3850, 3851, 3852, 3853, 3854, 3855, 3856, 3857, 3858, 3859, 3860, 3861, 3862, 3863, 3864, 3865, 3866, 3867, 3868, 3869, 3870, 3871, 3872, 3873, 3874, 3875, 3876, 3877, 3878, 3879, 3880, 3881, 3882, 3883, 3884, 3885, 3886, 3887, 3888, 3889, 3890, 3891, 3892, 3893, 3894, 3895, 3896, 3897, 3898, 3899, 3900, 3901, 3902, 3903, 3904, 3905, 3906, 3907, 3908, 3909, 3910, 3911, 3912, 3913, 3914, 3915, 3916, 3917, 3918, 3919, 3920, 3921, 3922, 3923, 3924, 3925, 3926, 3927, 3928, 3929, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 3937, 3938, 3939, 3940, 3941, 3942, 3943, 3944, 3945, 3946, 3947, 3948, 3949, 3950, 3951, 3952, 3953, 3954, 3955, 3956, 3957, 3958, 3959, 3960, 3961, 3962, 3963, 3964, 3965, 3966, 3967, 3968, 3969, 3970, 3971, 3972, 3973, 3974, 3975, 3976, 3977, 3978, 3979},
    {2048, 2049, 2050, 2051, 2052, 2053, 2054, 2057, 2058, 2061, 2062, 2065, 2066, 2067, 2068, 2070, 2071, 2072, 2075, 2076, 2079, 2080, 2081, 2082, 2085, 2087, 2089, 2091, 2092, 2094, 2096, 2097, 2099, 2100, 2101, 2102, 2105, 2106, 2107, 2110, 2111, 2113, 2114, 2116, 2117, 2119, 2120, 2125, 2126, 2127, 2128, 2129, 2130, 2132, 2133, 2134, 2135, 2137, 2138, 2141, 2143, 2144, 2145, 2146, 2147, 2148, 2149, 2151, 2152, 2153, 2155, 2156, 2157, 2158, 2159, 2160, 2161, 2162, 2163, 2166, 2167, 2168, 2169, 2170, 2171, 2173, 2174, 2175, 2176, 2177, 2178, 2180, 2181, 2182, 2183, 2185, 2186, 2187, 2189, 2191, 2192, 2193, 2194, 2195, 2196, 2197, 2198, 2199, 2202, 2203, 2204, 2205, 2207, 2208, 2209, 2211, 2214, 2215, 2217, 2219, 2221, 2222, 2223, 2224, 2225, 2226, 2227, 2228, 2230, 2231, 2232, 2234, 2235, 2236, 2238, 2239, 2241, 2242, 2243, 2245, 2246, 2247, 2248, 2249, 2250, 2251, 2252, 2254, 2255, 2256, 2259, 2260, 2262, 2263, 2264, 2265, 2267, 2268, 2269, 2270, 2272, 2273, 2274, 2277, 2279, 2280, 2281, 2282, 2283, 2284, 2285, 2286, 2288, 2289, 2291, 2292, 2293, 2294, 2295, 2296, 2297, 2298, 2301, 2302, 2303, 2304, 2305, 2310, 2312, 2313, 2314, 2315, 2316, 2317, 2318, 2319, 2320, 2321, 2322, 2325, 2327, 2328, 2330, 2332, 2333, 2335, 2336, 2337, 2338, 2341, 2342, 2343, 2344, 2345, 2346, 2347, 2349, 2350, 2351, 2355, 2357, 2358, 2360, 2361, 2362, 2364, 2365, 2366, 2367, 2368, 2371, 2372, 2373, 2375, 2377, 2378, 2379, 2380, 2382, 2383, 2384, 2385, 2387, 2388, 2389, 2390, 2391, 2393, 2394, 2396, 2397, 2398, 2399, 2400, 2401, 2402, 2403, 2405, 2406, 2411, 2412, 2413, 2415, 2416, 2417, 2419, 2420, 2421, 2422, 2424, 2425, 2426, 2427, 2429, 2431, 2432, 2434, 2435, 2436, 2437, 2438, 2439, 2440, 2441, 2443, 2444, 2445, 2447, 2448, 2449, 2450, 2451, 2452, 2453, 2454, 2455, 2456, 2457, 2458, 2459, 2461, 2463, 2465, 2466, 2468, 2470, 2471, 2472, 2473, 2474, 2475, 2476, 2478, 2479, 2480, 2481, 2483, 2486, 2487, 2488, 2490, 2491, 2493, 2494, 2496, 2497, 2498, 2501, 2502, 2503, 2505, 2508, 2509, 2510, 2511, 2512, 2513, 2514, 2515, 2516, 2517, 2518, 2519, 2522, 2523, 2525, 2527, 2528, 2529, 2530, 2531, 2533, 2534, 2535, 2537, 2538, 2540, 2541, 2543, 2544, 2545, 2547, 2548, 2555, 2557, 2558, 2562, 2565, 2566, 2567, 2568, 2569, 2570, 2571, 2572, 2576, 2577, 2580, 2581, 2582, 2583, 2584, 2585, 2587, 2588, 2589, 2592, 2594, 2595, 2596, 2597, 2598, 2599, 2603, 2605, 2608, 2609, 2610, 2612, 2614, 2616, 2617, 2618, 2620, 2621, 2622, 2626, 2627, 2628, 2629, 2632, 2633, 2634, 2635, 2636, 2637, 2639, 2640, 2641, 2642, 2643, 2644, 2645, 2647, 2648, 2649, 2650, 2651, 2652, 2653, 2656, 2657, 2658, 2659, 2660, 1465, 1577, 1718, 1912, 1913, 1914, 1915, 1916, 1919, 1920, 1921, 1922, 1923, 1924, 1926, 1927, 1928, 1930, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1939, 1940, 1941, 1942, 1944, 1945, 1947, 1948, 1949, 1950, 1951, 1952, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1964, 1965, 1968, 1969, 1970, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1980, 1981, 1982, 1987, 1988, 1990, 1991, 1992, 1994, 1995, 1996, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2031, 2032, 2034, 2035, 2036, 2038, 2039, 2041, 2042, 2044, 2047},
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 351, 364, 393, 399, 427, 441, 454, 464, 476, 501, 549, 564, 2704, 2740, 2814, 2838, 2885, 3003, 1171, 1193, 3290, 1297, 1387, 1486, 1549},
    {3073, 3078, 3080, 3084, 3085, 3090, 3092, 3095, 3096, 3098, 3104, 3108, 3114, 3119, 3123, 3129, 3130, 3135, 3137, 3139, 3142, 3143, 3149, 3153, 3158, 3163, 3166, 3174, 2663, 2664, 2666, 2669, 3184, 2673, 2676, 2679, 2680, 2681, 3193, 2683, 3195, 2689, 3202, 2694, 2695, 2698, 3213, 2705, 2706, 3220, 2710, 3226, 3228, 3229, 2717, 3231, 3232, 3234, 3237, 2726, 3238, 2729, 2731, 3243, 3250, 2738, 2741, 2743, 2745, 2746, 3257, 2749, 2750, 2754, 2755, 2756, 3269, 2757, 2761, 2763, 2766, 2769, 3284, 2773, 3286, 2777, 3293, 3294, 2782, 2786, 3300, 2789, 3301, 3303, 2794, 2798, 2800, 3315, 3317, 2806, 2807, 2809, 2810, 3324, 2815, 3333, 2821, 3339, 2827, 2835, 3348, 2837, 2839, 2850, 2851, 2852, 3362, 3363, 3364, 2854, 2862, 2864, 2866, 2867, 2872, 3384, 2873, 2874, 3388, 3387, 2877, 3390, 3392, 2880, 3394, 2884, 3399, 2888, 3400, 2887, 2890, 2891, 2896, 2897, 2899, 3413, 3414, 2905, 3417, 2906, 2908, 2907, 2909, 2910, 2911, 2912, 2913, 3426, 3428, 2915, 2917, 3430, 2916, 2919, 2920, 3433, 3435, 2924, 2925, 2927, 2928, 2929, 3436, 2931, 2940, 2943, 2944, 2945, 2948, 2953, 2956, 2960, 2963, 2966, 2967, 2969, 2974, 2977, 2978, 2981, 2985, 2986, 2987, 2988, 2993, 2994, 3000, 3004, 3010, 3014, 3022, 3023, 3025, 3026, 3029, 3030, 3033, 3036, 3040, 3042, 3044, 3045, 3047, 3052, 3054, 3056, 3059, 3064, 3065},
    {2560, 2561, 2563, 2564, 2055, 2056, 2059, 2060, 2573, 2574, 2575, 2063, 2064, 2578, 2579, 2069, 2073, 2586, 2074, 2549, 2077, 2590, 2591, 2078, 2593, 2083, 2084, 2086, 2600, 2088, 2601, 2602, 2090, 2604, 2606, 2093, 2095, 2607, 2098, 2611, 2613, 2103, 2104, 2615, 2619, 2108, 2109, 2623, 2112, 2624, 2625, 2115, 2630, 2631, 2118, 2121, 2122, 2123, 2124, 2638, 2131, 2646, 2136, 2139, 2140, 2654, 2142, 2655, 2150, 2154, 2164, 2165, 2172, 2179, 2184, 2188, 2190, 2200, 2201, 2206, 2210, 2212, 2213, 2216, 2218, 2220, 2229, 2233, 2237, 2240, 2244, 2552, 2253, 2257, 2258, 2261, 2266, 2271, 2275, 2276, 2278, 2287, 2290, 2299, 2300, 2306, 2307, 2308, 2309, 2311, 2323, 2324, 2326, 2329, 2331, 2334, 2339, 2340, 2348, 2352, 2353, 2354, 2356, 2359, 2363, 2369, 2370, 2374, 2376, 2381, 2386, 2392, 2395, 2404, 2407, 2408, 2409, 2410, 2414, 2418, 2423, 2428, 1917, 2430, 1918, 2433, 1925, 1929, 2442, 2446, 1938, 1943, 1946, 2460, 2462, 2464, 1953, 2467, 2469, 1962, 1963, 2477, 1966, 2482, 1971, 2484, 2485, 2489, 1979, 2492, 2495, 1983, 1984, 1985, 1986, 2499, 1989, 2500, 2504, 1993, 2506, 2507, 1997, 2005, 2520, 2521, 2524, 2526, 2532, 2021, 2020, 2536, 2539, 2542, 2030, 2033, 2546, 2037, 2550, 2551, 2040, 2553, 2554, 2043, 2556, 2045, 2046, 2559},
    {686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 858, 859, 860, 861, 863, 864, 866, 867, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895},
    {3980, 3981, 3982, 3983, 3984, 3985, 3986, 3987, 3988, 3989, 3990, 3991, 3992, 3993, 3994, 3995, 3996, 3997, 3998, 3999, 4000, 4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008, 4009, 4010, 4011, 4012, 4013, 4014, 4015, 4016, 4017, 4018, 4019, 4020, 4021, 4022, 4023, 4024, 4025, 4026, 4027, 4028, 4029, 4030, 4031, 4032, 4033, 4034, 4035, 4036, 4037, 4038},
    {901, 1798, 903, 1802, 1548, 1424, 1568, 1186, 1314, 1319, 938, 1452, 1837, 1333, 1468, 958, 1216, 1601, 963, 1733, 1097, 1232, 1747, 1363, 1493, 1880, 1629, 1118, 992, 1378, 1892, 1895, 1648, 1137, 1784, 1787, 1533},
    {640, 643, 647, 650, 658, 659, 661, 662, 665, 670, 675, 681, 576, 577, 578, 582, 583, 595, 599, 600, 615, 627, 628, 632, 635},{3008, 2982, 2699, 3309, 2670, 2959, 2767, 2703, 2834, 3283, 3311, 3314, 3382, 3318, 2879, 2972, 3325, 3423},{2817, 3074, 3147, 3055, 2774, 3127}};
    double mod = g.calculateModularity(comm);
    cout << mod << endl;
}

int main(){
    //distr_of_friends();
    //evolution_of_opinions();
    //degree_distr();
    Av_degree();
    //test();

   /* ofstream clusFile("Clustering_coefficient_WS_vs_beta.txt");
    for (int i = 0; i < 100; i++){
        double beta = double(i) / 100.;
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta);
        clusFile << beta << ' ' << g.overallClustering() << '\n'; 
    } 
    clusFile.close();*/

}



// maybe also include adjecency matrix

// QUESTION: does every class need a destructor? + can we assign edges in time slower than N^2?

// Standard variation --> welford's online algorithm (https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/) 

// TO DO: check echo chambers for SBM with p_add big, so looks more like ER --> needs further testing!

// for clustersizes drawn from power-law: http://antoineallard.github.io/graph_cpp_library/group__RandomNumberGenerators.html, http://antoineallard.github.io/graph_cpp_library/group__RandomNumberGenerators.html#ga273ce7be1b9abcf307fe726b7f7761ba