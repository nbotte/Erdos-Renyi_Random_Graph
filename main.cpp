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
 /*   int N = 1000;
    double p = 0.1;
    double p_bern = 0.1;
    double initOp0Frac = 0.5;*/
    //g.print();
    
   /* vector<double> fractionAt0(N);
    double oldFractionAt0;
    vector<double> fractionAt500(N);
    double oldFractionAt500;
    double binWidth = 0.1;
    int numberOfBins = 1/binWidth;
    vector<int> neighOp1HistAt0(numberOfBins);
    vector<int> neighOp1HistAt500(numberOfBins);
    // average over different network
    for (int n = 0; n < 10; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, initOp0Frac, 0);
        cout << "Graph " << n << endl;
        // for each network average over different simulation
        for (int s = 0; s < 10; s++){
            int opinion1 = 0;
            cout << "Simulation " << s << endl;
            // for each node count fraction of neighbors with opinion 1 at t = 0 (averaged)
            for (int i = 0; i < g.nodelist().size(); i++){
                if (g.nodelist()[i].neigh().size() != 0){
                    for (int index : g.nodelist()[i].neigh()){
                        if (g.nodelist()[index].opinion() == 1){
                            opinion1++;
                        }   
                    }
                }
                if (g.nodelist()[i].neigh().size() != 0){
                    oldFractionAt0 = double(opinion1)/g.nodelist()[i].neigh().size();
                }
                else{
                    oldFractionAt0 = 0.;
                }
                fractionAt0[i] += (oldFractionAt0/100.);
                opinion1 = 0;
            }

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
            g.resetInitOpinion();
        }
    }
    
    // make histogram for nodes with certain fraction of neighbors with opinion 1
    for (int i = 0; i < fractionAt500.size(); i++){
        // determine position of fraction in the histogram eg 0.111 lies in interval 0-0.1 so position is 0
        int index0 = int(fractionAt0[i]*10);
        int index500 = int(fractionAt500[i]*10);
        if (index0 == 10){
            index0--;
        }
        if (index500 == 10){
            index500--;
        }
        neighOp1HistAt0[index0] += 1;
        neighOp1HistAt500[index500] += 1;
    }

    vector<double> HistNorm(numberOfBins);
    // make normalized histogram at t = 500 by dividing by the one at t = 0
    for (int i = 0; i < numberOfBins; i++){
        double normVal;
        if (neighOp1HistAt0[i] != 0){
            normVal = double(neighOp1HistAt500[i])/double(neighOp1HistAt0[i]);
        }
        else{
            if (neighOp1HistAt500[i] == 0){
                normVal = 1.;
            }
            else{
                normVal = double(neighOp1HistAt500[i]);
            }
        }
        HistNorm[i] = normVal;
    }

    ofstream normfile("Normalized_hist_fraction_friends_opinion1_001_av.txt");
    for (int i = 0; i < HistNorm.size(); i++){
        normfile << HistNorm[i] << endl;
    }*/



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
    vector<int> clusterSizes(10);
    vector<double> edgeProbs(10);
    for (int i = 0; i < 10; i++){
        clusterSizes[i] = 100;
        edgeProbs[i] = 0.5;
    }

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
    /*ofstream opfile("Fraction_of_opinions_Clustered_01-001_50_50_no_stubb_paper8_active_01_av_good_init.txt");
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

    int N = 1000;
    int K = 20;
    double beta = 0.01;

    int count = 0;

    Watts_Strogatz_Network g = Watts_Strogatz_Network(N, K, beta); 
    double x = pow(1-beta, 3) * double(3*(K-2))/double((4*(K-1)));
    double x1 = double(K)/double(N-1);
    cout << g.overallClustering() << ' ' << x << endl;
   /* for (int i = 0; i < g.nodelist().size(); i++){
       // cout << g.nodelist()[i] << ": ";
        for (int index : g.nodelist()[i].neigh()){
            //cout << g.nodelist()[index] << ' ';
            count++;
        }
        //cout << endl;
    }

    int x = N * K/2;

    cout << count/2 << ' ' << x << endl;*/

    /*vector<int> degreeDistr(500);
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

// TO DO: add averaged histogram to shared file!!
// TO DO: check if Watts-Strogatz produces correct network (degree distribution etc --> see wiki: properties) + make graph of C(beta)/C(0) vs beta --> goes like (1-beta)^3?
