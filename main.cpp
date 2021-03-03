// Nina Botte

#include <cmath>
#include "ClusteredRandomER.h"
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
using namespace std;

// function that calculates the distribution of friends with the same opinion 1 for a particular network (according to paper 8)
void distr_of_friends(){
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
    double p_add = 0.0005;
    /*vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};*/
    vector<int> clusterSizes(10); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(10);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        edgeProbs.push_back(0.5);
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 100;
        edgeProbs[i] = 0.07;
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

    // average over different networks
    for (int n = 0; n < 10; n++){
        //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac);
        g.setNodeThreshold(0.8);
        //mod += (g.calculateModularity()/10.);
        
        cout << "Graph " << n << endl;

        // for each network average over different simulation
        for (int s = 0; s < 10; s++){
            cout << "Simulation " << s << endl;

            // reset the initial opinions to start a new simulation for the same network + make nodes stubborn
            g.resetInitOpinion(initOp0Frac);
            //g.makeRandomFractionStubborn(0.1);

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

    ofstream normfile("Hist_500_and_0_fraction_friends_opinion1_WS_PR_6-001_T=08.txt");
   // ofstream varfile("Hist_500_and_0_fraction_friends_opinion1_SBM_PR_01-0001_res=0_10x100_mean_var_11_bins.txt");
    ofstream xfile("Hist_500_and_0_fraction_friends_opinion1_WS_PR_6-001_T=08_xvalues.txt");
    ofstream echofile("Echo_chamber_WS_PR_6-001_T=08.txt");
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
    double p_add = 0.0005;
    /*vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};*/
    vector<int> clusterSizes(10); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(10);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        edgeProbs.push_back(0.5);
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 100;
        edgeProbs[i] = 0.06;
    }

    ofstream opfile("Fraction_of_opinions_WS_50_50_no_stubb_paper8_active_01_av_good_init_PR_6-001_T=08.txt");
    vector<double> mean0(500); // contains the average fraction of opinion 0 in the graph at each timestep
    vector<double> mean1(500); // contains the average fraction of opinion 1 in the graph at each timestep
    vector<double> variance0(500); // calculate variance of opinion 0 according to Welford's algorithm
    vector<double> variance1(500); // calculate variance of opinion 1 according to Welford's algorithm
    int count = 1;

    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 10; n++){
        //Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac);
        g.setNodeThreshold(0.8);
        cout << "Graph: " << n << endl;
        // for each network, run different simulations --> can this be implemented faster?
        for (int s = 0; s < 10; s++){
            // reset the initial opinions
            g.resetInitOpinion(initOp0Frac);
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
            count++;
        }
    }
    // for each time step: print the average opinion fraction over the different graphs to see the opinion evolution
    for (int i = 0; i < 500; i++){
        opfile << mean0[i] << ' ' << mean1[i] << ' ' << variance0[i]/double(count - 2) << ' ' << variance1[i]/double(count - 2) << endl;
    }
    opfile.close();
}

// function that calculates the degree distribution + the average degree distribution of a particular network
void degree_distr(){
    int N = 1000;
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

    vector<int> degreeDistr(500);
    Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
    for (int i = 0; i < g.nodelist().size(); i++){
        degreeDistr[g.nodelist()[i].neigh().size()] += 1;
    }
    ofstream degreeFile("Degree_distribution_Erdos-Renyi_01.txt");
    for (int i = 0; i < degreeDistr.size(); i++){
        degreeFile << i << ' ' << degreeDistr[i] << '\n';
    }

    degreeFile.close();
    vector<int> degreeDistrAv(500);
    for (int k = 0; k < 100; k++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
        cout << "Graph: " << k << endl;
        for (int i = 0; i < g.nodelist().size(); i++){
            degreeDistrAv[g.nodelist()[i].neigh().size()] += 1;
        }
    }


    ofstream degreeFileAv("Degree_distribution_Erdos-Renyi_01_av.txt");
    for (int i = 0; i < degreeDistrAv.size(); i++){
        degreeFileAv << i << ' ' << degreeDistrAv[i]/100 << '\n';
    }

    degreeFileAv.close();
}

// function that calculates the average degree of a graph
void Av_degree(){
    int N = 1000;
    int K = 10;
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
    double p_add = 0.007;
    /*vector<int> clusterSizes = {}; // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs = {};*/
    vector<int> clusterSizes(100); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(100);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        if (x < 100){
            edgeProbs.push_back(0.2);
        }
        else{
            edgeProbs.push_back(0.05);
        }        
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 10;
        edgeProbs[i] = 0.06;
    }
    int degree = 0;
    Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
    for (int i = 0; i < g.nodelist().size(); i++){
        degree += g.nodelist()[i].neigh().size();
    }
    double mod = g.calculateModularity();
    double Av_degree = double(degree)/double(N);
    cout << Av_degree << ' ' << mod << endl;
}

// function to do small tests
void test(){
    int N = 1000;
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
    double p_add = 0.002;
    vector<int> clusterSizes(50); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(50);
    /*double x;
    while (file >> x){
        clusterSizes.push_back(x);
        edgeProbs.push_back(0.5);
    }
    file.close();*/
    for (int i = 0; i < clusterSizes.size(); i++){
        clusterSizes[i] = 20;
        edgeProbs[i] = 0.25;
    }

    Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
    //Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac);
    cout << g.calculateModularity() << endl;
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
    /*g.setNodesActive(p_bern);
    g.setNodeThreshold(0.5);

    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index] << ' ' << g.nodelist()[index].threshold() << '\t';
        }
        cout << endl;
    }

    g.changeOpinions();
    g.changeOpinions();

    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index] << ' ' << g.nodelist()[index].threshold() << '\t';
        }
        cout << endl;
    }*/
}

int main(){
    distr_of_friends();
    evolution_of_opinions();
    //degree_distr();
    //Av_degree();
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