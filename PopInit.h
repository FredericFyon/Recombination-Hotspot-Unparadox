#ifndef POPINIT_H_INCLUDED
#define POPINIT_H_INCLUDED

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <chrono>

using namespace std;

vector<vector<vector<double> > > PopInit(int Nloc, int Npop, double pB_0, double pA_0, mt19937 mt)
{

    uniform_real_distribution<double> unifDouble(0.0,1.0);  //Continuous uniform drawing between 0 and 1
    uniform_int_distribution<int> unifInt(0,Npop-1);        //Discrete uniform drawing between 0 and Npop - 1

    vector<vector<vector<double> > > Pop;

    double p1 = Npop*pB_0*pB_0;                             //p1 = Number of B1 B1 individuals
    double p2 = p1 + Npop*pB_0*(1-pB_0);                    //p2 - p1 = Number of B1 B2 individuals
    double p3 = p2 + Npop*pB_0*(1-pB_0);                    //p3 - p2 = Number of B2 B1 individuals

    for(int i = 0; i < Npop; i++)
    {
        vector<vector<double> > IndivI;                     //Individual i

        //PRDM9 locus

        vector<double> LocPRDM9;                            //Locus A

        int Tirage_PRDM9_11 = unifInt(mt);                  //Drawing a integer between 0 and Npop - 1

        LocPRDM9.push_back(Tirage_PRDM9_11/(Npop/Nloc)+1);  //The first value of PRDM9 vector specifies which target locus
                                                            //the allele preferentially binds. In the following, we assume
                                                            //no preference such that this value is not used and does not matter

        double Tirage_PRDM9_12 = unifDouble(mt);            //Drawing a double between 0 and 1

        if(Tirage_PRDM9_12 < pA_0)                          //With probability pA_0, the first PRDM9 allele of the individual
        {                                                   //is A1, here coded as a 0
            LocPRDM9.push_back(0);
        }
        else                                                //With probability 1 - pA_0 it is A2, here coded as a 1
        {
            LocPRDM9.push_back(1);
        }

        int Tirage_PRDM9_21 = unifInt(mt);                  //Drawing a integer between 0 and Npop - 1

        LocPRDM9.push_back(Tirage_PRDM9_21/(Npop/Nloc)+1);  //Specifies to which target the second PRDM9 allele preferentially
                                                            //binds. Again, no use for that in the following

        double Tirage_PRDM9_22 = unifDouble(mt);            //Drawing a double between 0 and 1

        if(Tirage_PRDM9_22 < pA_0)                          //With probability pA_0, the second PRDM9 allele of the individual
        {                                                   //is A1, here coded as a 0
            LocPRDM9.push_back(0);
        }
        else                                                //With probability 1 - pA_0 it is A2, here coded as a 1
        {
            LocPRDM9.push_back(1);
        }

        IndivI.push_back(LocPRDM9);                         //The two PRDM9 alleles are added to the genome of individual i

        //Target loci

        for(int j = 0; j < Nloc; ++j)
        {
            vector<double> LocJ;                            //(j+1)th target locus

            int Tirage_LocJ = unifInt(mt);                  //Drawing a integer between 0 and Npop - 1

            if(Tirage_LocJ < p1)                            //The first p1 individuals are B1B1 at target locus j+1
            {
                LocJ.push_back(0);
                LocJ.push_back(0);
            }
            else if(Tirage_LocJ >= p1 && Tirage_LocJ < p2)  //Then, the following p2 - p1 individuals are B1B2 at target locus j+1
            {
                LocJ.push_back(0);
                LocJ.push_back(1);
            }
            else if(Tirage_LocJ >= p2 && Tirage_LocJ < p3) //Then, the following p3 - p2 individuals are B2B1 at target locus j+1
            {
                LocJ.push_back(1);
                LocJ.push_back(0);
            }
            else if(Tirage_LocJ >= p3 && Tirage_LocJ < Npop) //Finally, all remaining individuals are B2B2 at target locus j+1
            {
                LocJ.push_back(1);
                LocJ.push_back(1);
            }

            IndivI.push_back(LocJ);                         //Target locus j+1 is added to the genome of individual i
        }

        Pop.push_back(IndivI);                              //Individual i is added to the population
    }

    return Pop;                                             //The genome of all individuals in the population has been defined
}

#endif // POPINIT_H_INCLUDED

