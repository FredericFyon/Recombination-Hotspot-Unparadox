#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <chrono>

#include "PopInit.h"
#include "LifeCycle.h"

using namespace std;
using namespace chrono;

int main()
{
    //Parameter values:

    int Npop = 10000;                                   //population size
    int Nloc = 2;                                       //number of target loci
    int Ngen = 100000;                                  //total number of generations
    int Nit = 1;                                        //number of iterations
    int tBurn = 2000;                                   //number of generations of the burn-in phase
    int tSample = 10;                                   //sampling time: recording every 10 generations
    double b = 1.0;                                       //binding success rate in case of matching alleles in a binding attempt
    double r = 0.5;                                     //recombination rate between loci
    double c = 1;                                       //conversion rate
    double f = 0.4;                                     //fitness cost of the absence of Crossing-Over (CO)
    double mu_A = 0.0000001;                            //mutation rate on PRDM9 locus
    double mu_B = 0.0000001;                            //mutation rate on target loci
    double pB_0 = 0.95;                                 //initial frequency of B1
    double pA_0 = 0.99;                                 //initial frequency of A1

    for(int it = 0; it < Nit; ++it)
    {
        //Creation of the output file where the data is recorded
            //CO probability
        ostringstream fileResults;
        fileResults << "C:/Users/Fred/Dropbox/Hotspots/Mutation Model/RecHot_ManyTargets_FinalResults/probCOFinal_It" << it << "_c" << c << "_f" << f << "_Npop" << Npop << "_Nloc" << Nloc << "_muAsym" << mu << "_r" << r << ".dat";
        string varResults = fileResults.str();
        ofstream Results(varResults);
        Results<< setprecision(6);

            //A1 frequency
        ostringstream fileResultsPRDM9;
        fileResultsPRDM9 << "C:/Users/Fred/Dropbox/Hotspots/Mutation Model/RecHot_ManyTargets_FinalResults/freqPRDM9Final_It" << it << "_c" << c << "_f" << f << "_Npop" << Npop << "_Nloc" << Nloc << "_muAsym" << mu << "_r" << r << ".dat";
        string varResultsPRDM9 = fileResultsPRDM9.str();
        ofstream ResultsPRDM9(varResultsPRDM9);
        ResultsPRDM9 << setprecision(6);

            //B1 frequency
        ostringstream fileResultsLoc;
        fileResultsLoc << "C:/Users/Fred/Dropbox/Hotspots/Mutation Model/RecHot_ManyTargets_FinalResults/freqLocFinal_It" << it << "_c" << c << "_f" << f << "_Npop" << Npop << "_Nloc" << Nloc << "_muAsym" << mu << "_r" << r << ".dat";
        string varResultsLoc = fileResultsLoc.str();
        ofstream ResultsLoc(varResultsLoc);
        ResultsLoc << setprecision(6);

        vector<vector<vector<double> > > Pop;               //Pop is the total population

        for(int t = 0; t <= Ngen; ++t)
        {
            mt19937 mt(314475649*(it+1)*(t+1));             //Mersenne-Twister random number generator
            //cout << t << endl;
            if(t == 0)
            {
                Pop = PopInit(Nloc,Npop,pB_0,pA_0,mt);           //Initialization of the population
            }
            else if(t > 0 && t <= tBurn)
            {
               Pop = LifeCycle(Pop,r,mu_A,mu_B,b,0.,0.,Npop,Nloc,mt);    //Burn-in generations: f = 0 and c = 0
            }
            else if(t > tBurn)
            {
                Pop = LifeCycle(Pop,r,mu_A,mu_B,b,c,f,Npop,Nloc,mt);     //Normal generations
            }


            //Data recording every tSample generations
            if(t%tSample == 0)
            {
                cout << "Iteration : " << it + 1 << " t = " << t << ";" << endl;    //This to know where the simulation is at


                array<double, 10> freqLoc = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};        //Array of B1 frequencies. It can accomodate up to 10 loci

                ResultsLoc << t << " ";                                             //Creation of the data line corresponding to this generation

                for(int nloc = 0; nloc < Nloc; ++nloc)
                {
                    //Calculation of the CO probability in the population
                    double probCO(0.);
                    for(int n = 0; n < Npop; ++n)
                    {
                        int b11(0),b12(0),b21(0),b22(0);
                        if(Pop[n][0][1] == Pop[n][nloc+1][0])
                        {
                            b11 = b;
                        }

                        if(Pop[n][0][1] == Pop[n][nloc+1][1])
                        {
                            b12 = b;
                        }

                        if(Pop[n][0][3] == Pop[n][nloc+1][0])
                        {
                            b21 = b;
                        }

                        if(Pop[n][0][3] == Pop[n][nloc+1][1])
                        {
                            b22 = b;
                        }

                        probCO = probCO + (1./(4.*Nloc))*(b11+b12+b21+b22);


                        //Calculation of the B1 frequencies in the population
                        if(Pop[n][nloc+1][0] == 0)
                        {
                            freqLoc[nloc] = freqLoc[nloc] + 1;
                        }

                        if(Pop[n][nloc+1][1] == 0)
                        {
                            freqLoc[nloc] = freqLoc[nloc] + 1;
                        }
                    }

                    probCO = probCO/Npop;
                    freqLoc[nloc] = freqLoc[nloc]/(2.*Npop);

                    Results << probCO << " ";
                    ResultsLoc << freqLoc[nloc] << " ";

                }
                Results << endl;
                ResultsLoc<<endl;

                //Calculation of the frequency of A1 in the population

                double freqPRDM9 = 0.;

                for(int nind = 0; nind < Npop; ++nind)
                {
                    if(Pop[nind][0][1] == 0)
                    {
                        freqPRDM9 = freqPRDM9 + 1;
                    }

                    if(Pop[nind][0][3] == 0)
                    {
                        freqPRDM9 = freqPRDM9 + 1;
                    }
                }
                freqPRDM9 = freqPRDM9/(2*Npop);
                ResultsPRDM9 << t << " " << freqPRDM9 << endl;
            }
        }
    }
    return 0;
}
