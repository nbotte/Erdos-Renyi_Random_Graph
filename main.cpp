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

// look at static functions to make graph.utils that has makeER function 

int main(){
    // looks ok, but needs further testing
    // N = 1000 and p = 0.1 takes a really long time to run (>1h30min) --> unfortunately true...! --> but is better if you don't add edges to the edge list

    // needs more testing, does opinion dynamics works properly --> see paper 8, smart to make nodes active (+ update opinionlist) in constructor? 
   /* int N = 1000;
    int K = 6;
    double p = 0.1;
    double initOp0Frac = 0.5;
    double beta = 0.01; // should be small enough in order to deviate from random case
    double p_bern = 0.1;*/

    double p_bern = 0.1;
    int N = 1000;
    int p_add = 0.01;
    double initOp0Frac = 0.5;
    vector<int> clusterSizes(100); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(100);
    for (int i = 0; i < 100; i++){
        clusterSizes[i] = 10;
        edgeProbs[i] = 0.1;
    }

    vector<double> fractionAt0(N);
    vector<double> fractionAt500(N);
    double binWidth = 0.1;
    int numberOfBins = 1/binWidth;
    vector<int> neighOp1HistAt500(numberOfBins);
    vector<int> neighOp1HistAt0(numberOfBins);


    // average over different networks
    for (int n = 0; n < 10; n++){
        Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
        cout << "Graph " << n << endl;

        // for each network average over different simulation
        for (int s = 0; s < 10; s++){
            cout << "Simulation " << s << endl;

            // reset the initial opinions to start a new simulation for the same network
            g.resetInitOpinion(initOp0Frac);

            for (int i = 0; i < g.nodelist().size(); i++){
                int opinion1 = 0;
                if (g.nodelist()[i].neigh().size() != 0){
                    for (int index : g.nodelist()[i].neigh()){
                        if (g.nodelist()[index].opinion() == 1){
                            opinion1++;                    
                        }
                    }
                    fractionAt0[i] = double(opinion1)/double(g.nodelist()[i].neigh().size());
                    // determine position of fraction in the histogram eg 0.111 lies in interval 0-0.1 so position is 0
                    int index0 = int(fractionAt0[i]*10);
        
                    if (index0 == 10){
                        index0--;
                    }
                    neighOp1HistAt0[index0] += 1; 
                }
            }

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
                    // determine position of fraction in the histogram eg 0.111 lies in interval 0-0.1 so position is 0
                    int index500 = int(fractionAt500[i]*10);
                    if (index500 == 10){
                        index500--;
                    }
                    neighOp1HistAt500[index500] += 1;
                }
            }
        }
    } 

    //g.print();
    // for now: make histogram only for t=500, t=0 not necessary
  //  vector<double> fractionAt500(N);
    //double oldFractionAt500;
   /* for (int n = 0; n < 10; n++){
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac);
        cout << "Graph " << n << endl;
        // for each network average over different simulation
        for (int s = 0; s < 10; s++){
            int opinion1 = 0;
            cout << "Simulation " << s << endl;

            // let opinions evolve in time
            for (int i = 0; i < 500; i++){
                g.setNodesActive(p_bern);
                g.changeOpinions();
                g.deactivateNodes();
            }

            // count fraction of neighbors with opinion 1 at t = 500 + make histogram with number of nodes with a certain fraction of neighbors with opinion 1
            for (int i = 0; i < g.nodelist().size(); i++){
                if (g.nodelist()[i].neigh().size() != 0){
                    for (int index : g.nodelist()[i].neigh()){
                        if (g.nodelist()[index].opinion() == 1){
                            opinion1++;
                        }
                    }
                }
                if (g.nodelist()[i].neigh().size() != 0){
                    oldFractionAt500 = double(opinion1)/g.nodelist()[i].neigh().size();
                }
                else{
                    oldFractionAt500 = 0.;
                }
                fractionAt500[i] += (oldFractionAt500/100.);
                opinion1 = 0;
            }

            // reset the initial opinions to start a new simulation for the same network
            g.resetInitOpinion(initOp0Frac);
        }
    }*/

    ofstream normfile("Hist_500_and_0_fraction_friends_opinion1_Clustered_01-001_100-10_REC.txt");
    for (int i = 0; i < neighOp1HistAt500.size(); i++){
        double norm = double(neighOp1HistAt500[i]) / double(neighOp1HistAt0[i]);
        normfile << neighOp1HistAt500[i] << ' ' << neighOp1HistAt0[i] << ' ' << norm << endl;
    }
    normfile.close();



   /* for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (int j = 0; j < g.nodelist()[i].opinionlist().size(); j++){
            cout << g.nodelist()[i].opinionlist()[j] << ' ';
        }
        cout << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index];
        }
        cout << endl;
    }*/
  /*  g.setNodesActive(p_bern);
    g.changeOpinions();

    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (int j = 0; j < g.nodelist()[i].opinionlist().size(); j++){
            cout << g.nodelist()[i].opinionlist()[j] << ' ';
        }
        cout << ": ";
        for (int index : g._nodelist[i].neigh()){
            cout << g.nodelist()[index];
        }
        cout << endl;
    }*/


   /* vector<int> degreeDistr(500);
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
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0);
        cout << "Graph: " << k << endl;
        for (int i = 0; i < g.nodelist().size(); i++){
            degreeDistrAv[g.nodelist()[i].neigh().size()] += 1;
        }
    }


    ofstream degreeFileAv("Degree_distribution_Erdos-Renyi_01_av.txt");
    for (int i = 0; i < degreeDistrAv.size(); i++){
        degreeFileAv << i << ' ' << degreeDistrAv[i]/100 << '\n';
    }

    degreeFileAv.close();*/

    //ofstream op1file("Fraction_of_opinions_1_50_50_no_stubb_paper8_active_01_good_init.txt"); 
   /* Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0);  
    for (int t=0; t<1000; t++){
        g.setNodesActive(p_bern);
        cout << g.countOpinionFraction()[0] << ' ' << g.countOpinionFraction()[1] << endl;
        g.changeOpinions();
        g.deactivateNodes();
        //g.print();
        //cout << endl;
    }*/
  //  g.changeOpinions();
   // g.print();
   // op1file.close();*/

   /* ofstream opfile("Fraction_of_opinions_01_50_50_no_stubb_paper8_active_01_good_av_good_init.txt");
    vector<double> mean0(500); // contains the average fraction of opinion 0 in the graph at each timestep
    vector<double> mean1(500); // contains the average fraction of opinion 1 in the graph at each timestep
    vector<double> variance0(500); // calculate variance of opinion 0 according to Welford's algorithm
    vector<double> variance1(500); // calculate variance of opinion 1 according to Welford's algorithm
    int count = 1;
    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 10; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0); 
        cout << "Graph: " << n << endl;
        // for each network, run different simulations --> can this be implemented faster?
        for (int s = 0; s < 10; s++){
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
            // reset the initial opinions
            g.resetInitOpinion();
            count++;
        }
    }
    // for each time step: print the average opinion fraction over the different graphs to see the opinion evolution
    for (int i = 0; i < 500; i++){
        opfile << mean0[i] << ' ' << mean1[i] << ' ' << variance0[i]/double(count - 2) << ' ' << variance1[i]/double(count - 2) << endl;
    }
    opfile.close();*/

    /*double p_bern = 0.1;
    int N = 1000;
    int p_add = 0.001;
    double initOp0Frac = 0.5;
    vector<int> clusterSizes(50); // length of this vector determines the number of cluster and the elements determine the size of each cluster
    vector<double> edgeProbs(50);
    for (int i = 0; i < 50; i++){
        clusterSizes[i] = 20;
        edgeProbs[i] = 0.5;
    }*/

   /* Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add"); 
    g.print();
    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index] << ' ';
        }
        cout << endl;
    }*/

    /*Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add");
    cout << g.averageClustering() << endl;*/
   /* g.print();
    cout << g.countOpinionFraction()[0] << ' ' << g.countOpinionFraction()[1] << endl;
    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ' ' << g.nodelist()[i].neigh().size() << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index] << ' ';
        }
        cout << endl;
    }*/
   /* g.setNodesActive(p_bern);
    for (int t=0; t<1000; t++){
        cout << g.countOpinionFraction()[0] << ' ' << g.countOpinionFraction()[1] << endl;
        g.changeRandomOpinion();
        //g.print();
        //cout << endl;
    }*/

    // also calculate standard deviation here?
    /*ofstream opfile("Fraction_of_opinions_Clustered_Cluster4_05-0001_50_50_no_stubb_paper8_active_01_av_good_init.txt");
    vector<double> mean0(500); // contains the average fraction of opinion 0 in the graph at each timestep
    vector<double> mean1(500); // contains the average fraction of opinion 1 in the graph at each timestep
    vector<double> variance0(500); // calculate variance of opinion 0 according to Welford's algorithm
    vector<double> variance1(500); // calculate variance of opinion 1 according to Welford's algorithm
    int count = 1;

    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 10; n++){
        Clustered_Random_Network g = Clustered_Random_Network(N, clusterSizes, edgeProbs, p_add, "add"); 
        cout << "Graph: " << n << endl;
        // for each network, run different simulations --> can this be implemented faster?
        for (int s = 0; s < 10; s++){
            // reset the initial opinions
            g.resetInitOpinion(initOp0Frac);
            // for each network and each simulation: let the opinions evolve in time
            for (int t = 0; t < 500; t++){
                g.setNodesActive(p_bern);
                g.changeOpinions();

                double x0 = g.countOpinionFractionCluster(4)[0];
                double x1 = g.countOpinionFractionCluster(4)[1];

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
    opfile.close();*/
   /* int numberOfEdges = 0;
    for (int i = 0; i < g.nodelist().size(); i++){
      //  cout << g.nodelist()[i] << ": ";
     //   numberOfEdges = numberOfEdges + g.nodelist()[i].neigh().size();
        cout << g.localClustering(g._nodelist[i]) << ' ';
       // cout << endl;
    }
    cout << endl;
    //cout << numberOfEdges/2 << endl;
    cout << g.averageClustering() << endl;*/

   /* int N = 1000;
    int K = 20;
    double initOp0Frac = 0.5;
    double beta = 0.01; // should be small enough in order to deviate from random case
    double p_bern = 0.1;

    ofstream opfile("Fraction_of_opinions_WS_001-20_20_80_no_stubb_paper8_active_01_av_good_init.txt");
    vector<double> mean0(500); // contains the average fraction of opinion 0 in the graph at each timestep
    vector<double> mean1(500); // contains the average fraction of opinion 1 in the graph at each timestep
    vector<double> variance0(500); // calculate variance of opinion 0 according to Welford's algorithm
    vector<double> variance1(500); // calculate variance of opinion 1 according to Welford's algorithm
    int count = 1;

    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 10; n++){
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac); 
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
    opfile.close();*/

   /* ofstream clusFile("Clustering_coefficient_WS_vs_beta.txt");
    for (int i = 0; i < 100; i++){
        double beta = double(i) / 100.;
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta);
        clusFile << beta << ' ' << g.overallClustering() << '\n'; 
    } 
    clusFile.close();*/

    /*Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta, initOp0Frac); 

    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (int index : g.nodelist()[i].neigh()){
            cout << g.nodelist()[index] << ' ';
        }
        cout << endl;
    }*/

   /* vector<int> degreeDistr(500);
    for (int i = 0; i < g.nodelist().size(); i++){
        degreeDistr[g.nodelist()[i].neigh().size()] += 1;
    }
    ofstream degreeFile("Degree_distribution_Watts-Strogatz_05.txt");
    for (int i = 0; i < degreeDistr.size(); i++){
        degreeFile << i << ' ' << degreeDistr[i] << '\n';
    }

    degreeFile.close();
    vector<int> degreeDistrAv(500);
    for (int k = 0; k < 100; k++){
        Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta);
        cout << "Graph: " << k << endl;
        for (int i = 0; i < g.nodelist().size(); i++){
            degreeDistrAv[g.nodelist()[i].neigh().size()] += 1;
        }
    }


    ofstream degreeFileAv("Degree_distribution_Watts-Strogatz_05_av.txt");
    for (int i = 0; i < degreeDistrAv.size(); i++){
        degreeFileAv << i << ' ' << degreeDistrAv[i]/100 << '\n';
    }

    degreeFileAv.close();*/




};

// maybe also include adjecency matrix

// QUESTION: does every class need a destructor? + can we assign edges in time slower than N^2?

// Standard variation --> welford's online algorithm (https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/) 

// TO DO: check echo chambers for SBM with p_add big, so looks more like ER --> needs further testing!

// for clustersizes drawn from power-law: http://antoineallard.github.io/graph_cpp_library/group__RandomNumberGenerators.html, http://antoineallard.github.io/graph_cpp_library/group__RandomNumberGenerators.html#ga273ce7be1b9abcf307fe726b7f7761ba