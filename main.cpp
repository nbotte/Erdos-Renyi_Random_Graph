// Nina Botte

#include <cmath>
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
    // N = 1000 and p = 0.1 takes a really long time to run (>1h30min) --> not no more!
    int N = 1000;
    double p = 0.1;
    double p_bern = 1.; 
    //Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern);
    /*ofstream op1file("Fraction_of_opinions_1_80_20_one_op_stubb_50.txt");  
    for (int t=0; t<300; t++){
        cout << g.countOpinionFraction()[0] << ' ' << g.countOpinionFraction()[1] << endl;
        //g.print();
        //cout << endl;
        g.changeOpinions();
    }
    op1file.close();

    ofstream opfile("Fraction_of_opinions_1_70_30_no_stubb_bern_1_av.txt");
    vector<double> fractionsA(300);
    vector<double> fractionsB(300); 
    // loop over different networks to take averages of the fraction of opinions for each time step
    for (int n = 0; n < 100; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern); 
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
        opfile << fractionsA[i]/100. << ' ' << fractionsB[i]/100. << '\n';
    }
    opfile.close();*/

    // look at opinion fraction of neighbours of node i (see paper)
    // started implementing this, but not yet tested etc!
    int opinion0;
    vector<double> neighOpinions(300);
    for (int n = 0; n < 100; n++){
        Erdos_Renyi_Network g = Erdos_Renyi_Network(N, p, p_bern);
        for (int t = 0; t < 300; t++){
            for (Node* n : g.nodelist()[10].neigh()){
                // count number of neighbours with opinion 0
                if (n->opinion() == 0){
                    opinion0++;
                }
                }
            double oldOpinion0 =  neighOpinions[t];   
            neighOpinions[t] = oldOpinion0 + opinion0;
            opinion0 = 0;
            g.changeOpinions();
            }
        }
        
    }


};

// maybe also include adjecency matrix

// QUESTION: does every class need a destructor? + can we assign edges in time slower than N^2?


/* current status: opinion dynamics: takes often long time to run --> optimizations possible? ALSO: 50/50 case with no stubborn actors almost always leads to a stable 
situation of 54% of one opinion and 46% of the other, however one would expect a 50/50 situation
--> if you take more time steps and average over more simulations, this problem seems to disappear (still needs to be tested in some more depth)*/

// TO DO: look at opinion fraction of neighbours of node i (see paper) + start implementing clustered random networks