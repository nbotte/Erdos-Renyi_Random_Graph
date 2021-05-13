// Nina Botte -- Master thesis: Opinion dynamics on social networks with stubborn actors

#include <cmath>
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

#include <typeinfo>
using namespace std;

// function that calculates the distribution of friends with the same opinion 1 for a particular network (according to paper 8)
void distr_of_friends(){
    int N = 10680;
    int K = 4;
    double p = 0.01;
    double initOp0Frac = 0.5;
    double beta = 0.2; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    // part needed for clustered random network
    /*ifstream file;
    file.open("Communities_sizes_real_network_lastfm.txt", ios::in);
    if(!file){
        cout << "No such file";
    }*/
    double p_add = 0.00007;
    /*vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};
    vector<int> meanDegrees = {};*/
    vector<int> clusterSizes(60); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(60);
    vector<int> meanDegrees(60);
   /* double x;
    while (file >> x){
        int size = int(x);
        clusterSizes.push_back(x);
        if (x < 100){
            edgeProbs.push_back(0.05);
        }
        else{
            edgeProbs.push_back(0.01);
        }  
        //edgeProbs.push_back(0.025);  
        //meanDegrees.push_back(4);
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 178;
        edgeProbs[i] = 0.12;
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

    int xAt0op0 = 0;
    int xAt0op1 = 0;
    int xAt500op0 = 0;
    int xAt500op1 = 0;
    
    // vectors that contain the number of opinion 0 and opinion 1 echo chambers at t=0 and t=500 for each of the independent simulations --> used to calculate the standard deviation
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
    
    /*ifstream file;
    file.open("PGP.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    int z; // only needed if there are three columns in edgelist
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y >> z){
        // if nodes start from 1: x-1 and y-1; else if nodes start from 0: x and y
        edge[0] = x - 1;
        edge[1] = y - 1;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(N, edges);*/

    // average over different networks
    for (int n = 0; n < 10; n++){
        //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
        Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
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
            g.makeRandomFractionStubborn(0., 0.); // make all nodes resistant (change how resistant they are from 0, 1)
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
                    
                    // determine position of fraction in the histogram eg 0.111 lies in interval 0.1-0.2 so position is 1
                    int index0 = int(fractionAt0[i]*10);
                    if (index0 == 10){
                        index0--;
                    }
                    neighOp1HistAt0[index0] += 1;  
                }
            }

            xVecAt0op0.push_back(xAt0op0); // add number of opinion 0 echo chambers at t=0 for single simulation
            xVecAt0op1.push_back(xAt0op1); // add number of opinion 1 echo chambers at t=0 for single simulation

            xAt0op0 = 0; // reset number of opinion 0 echo chambers, so that you can calculate them for the next simulation
            xAt0op1 = 0; // reset number of opinion 1 echo chambers, so that you can calculate them for the next simulation

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

            xVecAt500op0.push_back(xAt500op0); // add number of opinion 0 echo chambers at t=500 for single simulation
            xVecAt500op1.push_back(xAt500op1); // add number of opinion 1 echo chambers at t=500 for single simulation
            
            xAt500op0 = 0; // reset number of opinion 0 echo chambers, so that you can calculate them for the next simulation
            xAt500op1 = 0; // reset number of opinion 1 echo chambers, so that you can calculate them for the next simulation

            count++; 
        }
    } 

    neighOp1HistAt0[0] = x_0_0;
    neighOp1HistAt0[neighOp1HistAt0.size() - 1] = x_0_1;
    neighOp1HistAt500[0] = x_500_0;
    neighOp1HistAt500[neighOp1HistAt500.size() - 1] = x_500_1;

    ofstream normfile("Hist_500_and_0_fraction_friends_opinion1_SBM-WS_REC_10680-4-012-000007_60x178_fracRes=0_stubb=0_av10x5.txt");
    ofstream xfile("Hist_500_and_0_fraction_friends_opinion1_SBM-WS_REC_10680-4-012-000007_60x178_fracRes=0_stubb=0_av10x5_xvalues.txt");
    ofstream echofile("Echo_chamber_SBM-WS_REC_10680-4-012-000007_60x178_fracRes=0_stubb=0_av10x5.txt");
    for (int i = 0; i < neighOp1HistAt500.size(); i++){
        double norm = double(neighOp1HistAt500[i]) / double(neighOp1HistAt0[i]);
        normfile << neighOp1HistAt500[i] << ' ' << neighOp1HistAt0[i] << ' ' << norm << endl;
    }

    for (int i = 0; i < 50; i++){
        xfile << xVecAt0op0[i] << ' ' << xVecAt0op1[i] << ' ' << xVecAt500op0[i] << ' ' << xVecAt500op1[i] << endl;
    }

    double echo_0 = double(x_500_0)/double(x_0_0);
    double echo_1 = double(x_500_1)/double(x_0_1);
    echofile << x_0_0 << ' ' << x_500_0 << ' ' << echo_0 << '\n' << x_0_1 << ' ' << x_500_1 << ' ' << echo_1 << endl;
    normfile.close();
    xfile.close();
    echofile.close();
    /*stubborn = stubborn/100;
    cout << stubborn << endl;*/

    //cout << mod << endl;
}

// function that calculates the evolution of opinions in a particular network
void evolution_of_opinions(){
    int N = 10680;
    int K = 4;
    double p = 0.01;
    double initOp0Frac = 0.5;
    double beta = 0.2; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    // part needed for clustered random network
    /*ifstream file;
    file.open("Community_sizes_powerlaw.txt", ios::in);
    if(!file){
        cout << "No such file";
    }*/
    double p_add = 0.00007;
    /*vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};*/
    vector<int> clusterSizes(60); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(60);
    vector<int> meanDegrees(60);
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
        clusterSizes[i] = 178;
        edgeProbs[i] = 0.12;
        meanDegrees[i] = 4;
    }
  
    ofstream opfile("Fraction_of_opinions_active_01_av_good_init_SBM-WS_REC_10680-4-012-000007_60x178_fracRes=0_stubb=0_av10x5.txt");
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

    /*ifstream file;
    file.open("PGP.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    int z; // only needed if there are three columns in edgelist
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y >> z){
        // if nodes start from 1: x-1 and y-1; else if nodes start from 0: x and y
        edge[0] = x - 1;
        edge[1] = y - 1;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(N, edges);*/

    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 10; n++){
       // Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
        Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
        //Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac, 0);
        //Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
       // g.setNodeThreshold(0.);
        cout << "Graph: " << n << endl;
        // for each network, run different simulations --> can this be implemented faster?
        for (int s = 0; s < 5; s++){
            // reset the initial opinions
            g.resetInitOpinion(initOp0Frac);
          //  g.setNodeThreshold(0.);
            g.makeRandomFractionStubborn(0., 0.); // make all nodes resistant (change how resistant they are from 0, 1)

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
    int N = 10680;
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
    file.open("PGP.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    int z; // only needed if there are three columns in edgelist
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y >> z){
        // if nodes start from 1: x-1 and y-1; else if nodes start from 0: x and y
        edge[0] = x - 1;
        edge[1] = y - 1;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(N, edges);

    vector<int> degreeDistr(150);
   // Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i].neigh().size() << endl;
        degreeDistr[g.nodelist()[i].neigh().size()] += 1;
    }
    ofstream degreeFile("Degree_distribution_real_network_PGP.txt");
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
    int N = 10680;
    int K = 4;
    double p = 0.01;
    double initOp0Frac = 0.5;
    double beta = 0.2; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    // part needed for clustered random network
    /*ifstream file;
    file.open("Communities_sizes_real_network_lastfm.txt", ios::in);
    if(!file){
        cout << "No such file";
    }*/
    double p_add = 0.00005;
   /* vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};
    vector<int> meanDegrees = {};*/
    vector<int> clusterSizes(60); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(60);
    vector<int> meanDegrees(60);

/*    double x;
    while (file >> x){
        int size = int(x);
        clusterSizes.push_back(x);
        if (x < 100){
            edgeProbs.push_back(0.025);
        }
        else{
            edgeProbs.push_back(0.025);
        }  
        edgeProbs.push_back(0.025);  
        meanDegrees.push_back(4);    
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 178;
        edgeProbs[i] = 0.02;
        meanDegrees[i] = 4;
    }
   /* ifstream file;
    file.open("PGP.txt", ios::in);
    if(!file){
        cout << "No such file";
    }

    int x;
    int y;
    int z; // only needed if there are three columns in edgelist
    vector<int> edge(2);
    vector<vector<int>> edges = {};
    while (file >> x >> y >> z){
        // if nodes start from 1: x-1 and y-1; else if nodes start from 0: x and y
        edge[0] = x-1;
        edge[1] = y-1;
        edges.push_back(edge);
    }
    file.close();

    Real_World_Network g = Real_World_Network(N, edges);*/

    //int degree = 0;
    //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
    //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
    //Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac, 0);
    /*for (int i = 0; i < g.nodelist().size(); i++){
        degree += g.nodelist()[i].neigh().size();
    }
    double mod = g.calculateModularity();
    //double modTest = g.calculateModularityTest(g.clusters());
    double Av_degree = double(degree)/double(N);
    double clus = g.averageClustering();*/
    ofstream clusFile("avDeg_SBM.txt");
    for (int i = 0; i < 10; i++){
        cout << i << endl;
       // Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac, 0);
        Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
       // Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
        //Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
        //double mod = g.calculateModularity();
        //double clus = g.averageClustering();
        int degree = 0;
        for (int i = 0; i < g.nodelist().size(); i++){
            degree += g.nodelist()[i].neigh().size();
        }
        double Av_degree = double(degree)/double(N);
        clusFile << Av_degree << '\n';      
    }
    clusFile.close();

    //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, meanDegrees, p_add, "add");
    /*int edges_tot = g.numberOfEdges();
    cout << Av_degree << ' ' << mod << ' ' << clus << ' ' << edges_tot << endl;*/
}

// function to do small tests
void test(){
    int N = 8003;
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
}

int main(){
    distr_of_friends();
    evolution_of_opinions();
    //degree_distr();
    //Av_degree();
   // test();

   /* ofstream clusFile("Clustering_coefficient_WS_vs_beta.txt");
    for (int i = 0; i < 100; i++){
        double beta = double(i) / 100.;
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta);
        clusFile << beta << ' ' << g.overallClustering() << '\n'; 
    } 
    clusFile.close();*/

}

// Standard variation --> welford's online algorithm (https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/) 

