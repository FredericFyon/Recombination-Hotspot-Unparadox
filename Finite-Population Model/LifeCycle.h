#ifndef LIFECYCLE_H_INCLUDED
#define LIFECYCLE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <random>
#include <vector>
#include <time.h>
#include <chrono>

#include "Meiosis.h"

using namespace std;
using namespace chrono;

vector<vector<vector<double> > > LifeCycle(vector<vector<vector<double> > > Pop, double R, double mu_A, double mu_B, double b, double c, double f, int Npop, int Nloc, mt19937 mt)
{
    vector<vector<vector<double> > > PopPr;                                                //Population at the next generation
    uniform_real_distribution<double> unifDouble(0.0,1.0);                                 //Continuous uniform drawing between 0 and 1
    uniform_int_distribution<int> unifInt(0,Npop-1);                                       //Discrete uniform drawing between 0 and Npop - 1

    vector<vector<double> > Chrom_Parent1,Chrom_Parent2;                                    //Transmitted chromosomes from parents 1 and 2

    for(int i = 0; i < Npop; ++i)
    {
        vector<vector<double> > Genome_OffspringI;                                          //Genome of offspring i of the next-generation population

        //Finding parents to offspring i

        bool Accepted1 = false;
        int Parent1, Parent2;

        while(Accepted1 == false)                                                           //This following loop is iterated until a parent 1 is
        {                                                                                   //accepted

            Parent1 = unifInt(mt);                                                          //Drawing parent 1

            double Selec = unifDouble(mt);                                                  //Drawing a double between 0 and 1

            vector<vector<vector<double> > > MeiosisPar1 = Meiosis(Pop[Parent1],R,mu_A,mu_B,b,c,f,Nloc,mt);    //Running the meiosis of parent 1

            double Fit = MeiosisPar1[0][0][0];                                              //The fitness of parent 1 is retrieved

            if(Selec < Fit)                                                                 //With probability equal to the fitness, parent 1 is
            {                                                                               //accepted
                Chrom_Parent1 = MeiosisPar1[1];                                             //In this case, we retrieve the chromosome parent 1 transmits
                Accepted1 = true;                                                           //to offspring i
            }
        }

        //The same algorithm is repeated to draw and accept parent 2 based on its fitness
        bool Accepted2 = false;
        while(Accepted2 == false)
        {
            Parent2 = unifInt(mt);

            double Selec2 = unifDouble(mt);

            vector<vector<vector<double> > > MeiosisPar2 = Meiosis(Pop[Parent2],R,mu_A,mu_B,b,c,f,Nloc,mt);
            double Fit2 = MeiosisPar2[0][0][0];

            if(Selec2 < Fit2)
            {
                Chrom_Parent2 = MeiosisPar2[1];                                             //The chromosome parent 2 transmits to offspring i is determined
                Accepted2 = true;
            }
        }

        //PRDM9 locus of offspring i
        vector<double> OffspringI_PRDM9;
        OffspringI_PRDM9.push_back(Chrom_Parent1[0][0]);
        OffspringI_PRDM9.push_back(Chrom_Parent1[0][1]);
        OffspringI_PRDM9.push_back(Chrom_Parent2[0][0]);
        OffspringI_PRDM9.push_back(Chrom_Parent2[0][1]);

        Genome_OffspringI.push_back(OffspringI_PRDM9);                          //PRDM9 locus is added to the genome of offspring i

        //Target locus j+1 of offspring i
        for(int j = 0; j < Nloc; ++j)
        {
            vector<double> OffspringI_LocJ;

            OffspringI_LocJ.push_back(Chrom_Parent1[j+1][0]);
            OffspringI_LocJ.push_back(Chrom_Parent2[j+1][0]);

            Genome_OffspringI.push_back(OffspringI_LocJ);                       //Target locus j+1 is added to offspring i
        }

        PopPr.push_back(Genome_OffspringI);                                     //Offspring i is added to the next-generation population
    }

    return PopPr;
}


#endif // LIFECYCLE_H_INCLUDED
