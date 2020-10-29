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
    // This example seem to work as expected
   /* Node n = Node(1);
    Node m = Node(2);
    Node k = Node(0);
    Node l = Node(3);

    Edge e1 = Edge(&n, &m);
    Edge e2 = Edge(&n, &k);
    Edge e3 = Edge(&n, &l);
    Edge e4 = Edge(&m, &k);
    Edge e5 = Edge(&l, &k);

    n.addNeigh(m.index());
    m.addNeigh(n.index());
    n.addNeigh(k.index());
    k.addNeigh(n.index());
    n.addNeigh(l.index());
    l.addNeigh(n.index());
    m.addNeigh(k.index());
    k.addNeigh(m.index());
    l.addNeigh(k.index());
    k.addNeigh(l.index());

    n.removeNeigh(m.index());

    vector<int> neigh = n.neigh();
    for (int i = 0; i < neigh.size(); i++){
        cout << neigh[i] << endl;   
    }*/

    // looks ok, but needs further testing
    // N = 1000 and p = 0.1 takes a really long time to run (>1h30min) --> unfortunately true...! --> but is better if you don't add edges to the edge list

    // needs more testing, does opinion dynamics works properly??
    int N = 1000;
    double p = 0.1;
    double p_bern = 1.;
    /*Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0);
    g.print();*/
    /*for (int i = 0; i < g.nodelist().size(); i++){
        cout << g.nodelist()[i] << ": ";
        for (Node* n : g.nodelist()[i].neigh()){
            cout << *n;
        }
        cout << endl;
    }*/
   
    /*ofstream op1file("Fraction_of_opinions_1_80_20_one_op_stubb_50.txt");  
    for (int t=0; t<300; t++){
        cout << g.countOpinionFraction()[0] << ' ' << g.countOpinionFraction()[1] << endl;
        //g.print();
        //cout << endl;
        g.changeOpinions();
    }*/
  //  g.changeOpinions();
   // g.print();
    //op1file.close();

  // ofstream opfile("Fraction_of_opinions_1_50_50_no_stubb_bern_025_av.txt");
    vector<double> fractionsA(300);
    vector<double> fractionsB(300); 
    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 100; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0); 
        // for each network: let the opinions evolve in time
        for (int t = 0; t < 300; t++){
            double oldFractionA = fractionsA[t];
            double oldFractionB = fractionsB[t];
            fractionsA[t] = oldFractionA + g.countOpinionFraction()[0];
            fractionsB[t] = oldFractionB + g.countOpinionFraction()[1];
            g.changeOpinions();
        }
    }
    // for each time step: print the average opinion fraction over the different graphs to see the opinion evolution
    for (int i = 0; i < 300; i++){
        cout << fractionsA[i]/100. << ' ' << fractionsB[i]/100. << endl;
    }
  //  opfile.close();

    // NOTE: this takes long time to run!
    // TO DO: test this further
    // count for each node the average fraction of neighbours with opinion 0 at t=299 and at t=0
    /*int opinion0Begin;
    int opinion0End;
    vector<double> numberOfNeigh0At0(N);
    vector<double> numberOfNeigh0At299(N);
    double oldFractionBegin;
    double oldFractionEnd;
    for (int n = 0; n < 100; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern, 0);
        // count for t = 0
        for (int i = 0; i < N; i++){
            for (Node* n : g.nodelist()[i].neigh()){
                // count number of neighbours with opinion 0
                if (n->opinion() == 0){
                    opinion0Begin++;
                }
            }
            oldFractionBegin = numberOfNeigh0At0[i];
            numberOfNeigh0At0[i] = oldFractionBegin + (double(opinion0Begin)/double(g.nodelist()[i].neigh().size()));
            opinion0Begin = 0;
        }  
        for (int t = 0; t < 300; t++){
            g.changeOpinions();
        }   
        // count for t = 299 
        for (int i = 0; i < N; i++){
            for (Node* n : g.nodelist()[i].neigh()){
                // count number of neighbours with opinion 0
                if (n->opinion() == 0){
                    opinion0End++;
                }
            }
            oldFractionEnd = numberOfNeigh0At299[i];
            numberOfNeigh0At299[i] = oldFractionEnd + (double(opinion0End)/double(g.nodelist()[i].neigh().size()));
            opinion0End = 0;
        }    
    }

    ofstream His0File("histogramAt0_20/80.txt");
    ofstream His299File("histogramAt299_20/80.txt");
    ofstream HisNormFile("normalizedHist_20/80.txt");

    // make a histogram out of the vectors numberOfNeigh
    int numberOfBins = 11;
    vector<int> histogramAt0(numberOfBins); // declare an empty histogram with size 10 (for t = 0)
    vector<int> histogramAt299(numberOfBins); // declare an empty histogram with size 10 (for t = 299)
    for (int i = 0; i < N; i++){
        double r0 = numberOfNeigh0At0[i]/100;
        double r299 = numberOfNeigh0At299[i]/100;
        if (r0 < 0.05){
            histogramAt0[0] += 1;
        }
        else if (0.05 <= r0 && r0 < 0.15){
            histogramAt0[1] += 1;
        }
        else if (0.15 <= r0 && r0 < 0.25){
            histogramAt0[2] += 1;
        }
        else if (0.25 <= r0 && r0 < 0.35){
            histogramAt0[3] += 1;
        }
        else if (0.35 <= r0 && r0 < 0.45){
            histogramAt0[4] += 1;
        }
        else if (0.45 <= r0 && r0 < 0.55){
            histogramAt0[5] += 1;
        }
        else if (0.55 <= r0 && r0 < 0.65){
            histogramAt0[6] += 1;
        }
        else if (0.65 <= r0 && r0 < 0.75){
            histogramAt0[7] += 1;
        }
        else if (0.75 <= r0 && r0 < 0.85){
            histogramAt0[8] += 1;
        }
        else if (0.85 <= r0 && r0 < 0.95){
            histogramAt0[9] += 1;
        }
        else if (0.95 <= r0){
            histogramAt0[10] += 1;
        }
        if (r299 < 0.05){
            histogramAt299[0] += 1;
        }
        else if (0.05 <= r299 && r299 < 0.15){
            histogramAt299[1] += 1;
        }
        else if (0.15 <= r299 && r299 < 0.25){
            histogramAt299[2] += 1;
        }
        else if (0.25 <= r299 && r299 < 0.35){
            histogramAt299[3] += 1;
        }
        else if (0.35 <= r299 && r299 < 0.45){
            histogramAt299[4] += 1;
        }
        else if (0.45 <= r299 && r299 < 0.55){
            histogramAt299[5] += 1;
        }
        else if (0.55 <= r299 && r299 < 0.65){
            histogramAt299[6] += 1;
        }
        else if (0.65 <= r299 && r299 < 0.75){
            histogramAt299[7] += 1;
        }
        else if (0.75 <= r299 && r299 < 0.85){
            histogramAt299[8] += 1;
        }
        else if (0.85 <= r299 && r299 < 0.95){
            histogramAt299[9] += 1;
        }
        else if (0.95 <= r299){
            histogramAt299[10] += 1;
        }
    }
    vector<double> normalizedHist(numberOfBins);
    for (int i = 0; i < 11; i++){
        if (histogramAt0[i] == 0){
            if (histogramAt299[i] == 0){
                normalizedHist[i] = 1.;
            }
            else{
                normalizedHist[i] = 1. + double(histogramAt299[i])/double(N);
            }
        }
        else{
            if (histogramAt299[i] == 0){
                normalizedHist[i] = 1. - double(histogramAt0[i])/double(N);
            }
            else{
                normalizedHist[i] = double(histogramAt299[i])/double(histogramAt0[i]);
                cout << normalizedHist[i] << ' ';
            }
        }
    }

    for (int i = 0; i < 11; i++){
        His0File << histogramAt0[i] << '\n';
        His299File << histogramAt299[i] << '\n';
        HisNormFile << normalizedHist[i] << '\n';
    }
    His0File.close();
    His299File.close();
    HisNormFile.close();*/

    // still needs copy constructor for nodes
 /*   Clustered_Random_Network g = Clustered_Random_Network(1.);
    g.print();
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

// TO DO: look if you can omit the edge list, because adding the edges takes a long time when working with shared_ptr
