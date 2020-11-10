// Nina Botte

#include <cmath>
#include "ClusteredRandomER.h"
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
using namespace std;

// look at static functions to make graph.utils that has makeER function 

int main(){
    // looks ok, but needs further testing
    // N = 1000 and p = 0.1 takes a really long time to run (>1h30min) --> unfortunately true...! --> but is better if you don't add edges to the edge list

    // needs more testing, does opinion dynamics works properly --> see paper 8, smart to make nodes active (+ update opinionlist) in constructor? 
    int N = 1000;
    double p = 0.01;
    double p_bern = 1.;
    Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0);
   // g.print();

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
    }
    g.setNodesActive(p_bern);
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


    vector<int> degreeDistr(500);
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

    degreeFileAv.close();

   /* ofstream op1file("Fraction_of_opinions_1_50_50_no_stubb_paper8_active_1_one_node.txt"); 
    Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0);  
    g.setNodesActive(p_bern);
    for (int t=0; t<10000; t++){
        op1file << g.countOpinionFraction()[0] << ' ' << g.countOpinionFraction()[1] << endl;
        g.changeRandomOpinion();
        //g.print();
        //cout << endl;
    }
  //  g.changeOpinions();
   // g.print();
    op1file.close();

    ofstream opfile("Fraction_of_opinions_1_50_50_no_stubb_paper8_active_1_one_node_av.txt");
    vector<double> fractionsA(10000);
    vector<double> fractionsB(10000); 
    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 100; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0); 
        g.setNodesActive(p_bern);
        cout << "Graph: " << n << endl;
        // for each network: let the opinions evolve in time
        for (int t = 0; t < 10000; t++){
            g.changeRandomOpinion();
            double oldFractionA = fractionsA[t];
            double oldFractionB = fractionsB[t];
            fractionsA[t] = oldFractionA + g.countOpinionFraction()[0];
            fractionsB[t] = oldFractionB + g.countOpinionFraction()[1];
        }
    }
    // for each time step: print the average opinion fraction over the different graphs to see the opinion evolution
    for (int i = 0; i < 10000; i++){
        opfile << fractionsA[i]/100. << ' ' << fractionsB[i]/100. << endl;
    }
    opfile.close();*/


   /* Clustered_Random_Network g = Clustered_Random_Network(1.);
    
    for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (shared_ptr<Node> n : g.nodelist()[i].neigh()){
            cout << *n;
        }
        cout << endl;
    }*/
};

// maybe also include adjecency matrix

// QUESTION: does every class need a destructor? + can we assign edges in time slower than N^2?

// Question: opinion dynamics in paper 8 has a less deterministic majority model: if 1/3 of opinions is 0, you change to opinion 0 with a probability of 1/3 --> do I need to change to this model?

// TO DO: ckeck opinion dynamics really carefull, correctly implemented? Do some more testing --> also maybe send mail!
// Don't get a stable configuration (no equilibrium)?? --> what if you start with perfect initial conditions? --> if you run for 5000 timesteps, you keep hovering around the initial distribution (good?)
