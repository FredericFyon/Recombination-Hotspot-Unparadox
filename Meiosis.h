#ifndef MEIOSIS_H_INCLUDED
#define MEIOSIS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <chrono>

using namespace std;
using namespace chrono;

vector<vector<vector<double> > >Meiosis(vector<vector<double> > Indiv, double R, double mu_A, double mu_B, double b, double c, double f, int Nloc, mt19937 mt, ofstream &Stats)
{
    vector<vector<vector<double> > > TransmChromAndFit;                     //A vector made of the fitness of the individual, and of the chromosome it
                                                                            //will transmit if parent of an offspring

    vector<vector<double> > TransmittedChromosome;                          //Transmitted chromosome to offspring

    uniform_real_distribution<double> unifDouble(0.,1.);                    //Continuous uniform drawing between 0 and 1
    uniform_int_distribution<int> unifInt(0,Nloc-1);                        //Discrete uniform drawing between 0 and Nloc - 1
    uniform_int_distribution<int> unifInt2(0,2*Nloc-1);                     //Discrete uniform drawing between 0 and 2*Nloc - 1

    vector<vector<double> > W;                                              //Fitness of the individual. Must be a vector of vector (2D table) to fit within
                                                                            //TransmChromAndFit (which is a vector of vector of vector(3D table))

    vector<int> RecordChrom1;                                               //Record which chromosome is type I meiotic product
    vector<int> RecordChrom2;                                               //Record which chromosome is type II meiotic product

    vector<vector<double> > IndivBuffer = Indiv;                            //Genome of the individual after recombination

    RecordChrom1.push_back(0);                                              //Type I starts at chromosome 1 (0)
    RecordChrom2.push_back(1);                                              //Type II starts at chromosome 2(1)

    for(int i = 1; i < Nloc+1; ++i)
    {
        double sampleR = unifDouble(mt);                                    //Drawing a double between 0 and 1

        if(sampleR < R)                                                     //With probability R there is a recombination event between locus i and i+1
        {
            IndivBuffer[i][0] = Indiv[i][abs(RecordChrom1[i-1]-1)];         //If recombination, allele of type I meiotic product of locus i is taken from the
                                                                            //other chromosome that was taken allele allele of locus i - 1
                                                                            //For example, if allele of type I at locus 3 was taken from
                                                                            //first chromosome, then allele at locus 4 is taken from second chromosome

            IndivBuffer[i][1] = Indiv[i][abs(RecordChrom2[i-1]-1)];         //The same thing is done for type II meiotic product

            RecordChrom1.push_back(abs(RecordChrom1[i-1]-1));               //Record that now type I is taken from the other chromosome
            RecordChrom2.push_back(abs(RecordChrom2[i-1]-1));               //Same for type II
        }
        else                                                                //If there is no recombination, allele at locus i is taken from the same chromosome
        {                                                                   //than allele at locus i - 1
            IndivBuffer[i][0] = Indiv[i][RecordChrom1[i-1]];
            IndivBuffer[i][1] = Indiv[i][RecordChrom2[i-1]];
            RecordChrom1.push_back(RecordChrom1[i-1]);                      //And again this is recorded
            RecordChrom2.push_back(RecordChrom2[i-1]);
        }
    }

    double attempt = unifDouble(mt);                                        //Drawing a double between 0 and 1
    double success = unifDouble(mt);                                        //Drawing a double between 0 and 1
    int tirageLocTarget = unifInt(mt);                                      //Drawing an integer between 0 and Nloc - 1: this is the target locus that
                                                                            //a PRDM9 will attempt to bind

    if(trial < 1./4.)                                                       //With probability 1/4 a given PRDM9 allele attempts to bind a given target allele
                                                                            //at the drawn target locus
    {
        if((IndivBuffer[0][1] == IndivBuffer[tirageLocTarget+1][0]) && succes <= b)     //If the PRDM9 allele 1 and target allele 1 match, and with probability b
        {                                                                               //the binding attempt is successful
            vector<double> WW;                                              //A fitness vector to be added to the 2D fitness table
            WW.push_back(1.);                                               //There has been a CO: fitness of the individual is 1
            W.push_back(WW);                                                //Fitness is added to the 2D fitness table

            int tirageConv = unifDouble(mt);                                //Drawing a double between 0 and 1

            if(tirageConv < c)                                              //With probability c there is a conversion following the CO
            {
                IndivBuffer[tirageLocTarget+1][0] = IndivBuffer[tirageLocTarget+1][1];  //In that case, the bound target is converted into its homolog
            }
        }
        else                                                                //If the binding attempt is not successful
        {
            vector<double> WW;
            WW.push_back(1. - f);                                           //Fitness of the individual is 1 - f
            W.push_back(WW);
        }
    }

    //The same code above to determine binding attempt, success, fitness and conversion is repeated for the 4 possible binding attempts between a PRDM9
    //allele and a target allele at the drawn target locus

    else if(trial >= 1./4. && trial < 1./2.)
    {
        if((IndivBuffer[0][1] == IndivBuffer[tirageLocTarget+1][1]) && succes <= b)     //PRDM9 allele 1 and target allele 2
        {
            vector<double> WW;
            WW.push_back(1.);
            W.push_back(WW);
            int tirageConv = unifDouble(mt);
            if(tirageConv < c)
            {
                IndivBuffer[tirageLocTarget+1][1] = IndivBuffer[tirageLocTarget+1][0];
            }
        }
        else
        {
            vector<double> WW;
            WW.push_back(1. - f);
            W.push_back(WW);
        }
    }
    else if(trial >= 1./2. && trial < 3./4.)
    {
        if((IndivBuffer[0][3] == IndivBuffer[tirageLocTarget+1][0]) && succes <= b)     //PRDM9 allele 2 and target allele 1
        {
            vector<double> WW;
            WW.push_back(1.);
            W.push_back(WW);
            int tirageConv = unifDouble(mt);
            if(tirageConv < c)
            {
                IndivBuffer[tirageLocTarget+1][0] = IndivBuffer[tirageLocTarget+1][1];
            }
        }
        else
        {
            vector<double> WW;
            WW.push_back(1. - f);
            W.push_back(WW);
        }
    }
    else
    {
        if((IndivBuffer[0][3] == IndivBuffer[tirageLocTarget+1][1]) && succes <= b)     //PRDM9 allele 2 and target allele 2
        {
            vector<double> WW;
            WW.push_back(1.);
            W.push_back(WW);
            int tirageConv = unifDouble(mt);
            if(tirageConv < c)
            {
                IndivBuffer[tirageLocTarget+1][1] = IndivBuffer[tirageLocTarget+1][0];
            }
        }
        else
        {
            vector<double> WW;
            WW.push_back(1. - f);
            W.push_back(WW);
        }
    }

    TransmChromAndFit.push_back(W);                                         //The fitness 2D table is the first element to be added to the TransmChromAndFit
                                                                            //3D table

    //Modelling segregation
    double Segreg = unifDouble(mt);                                         //Drawing a double between 0 and 1

    if(Segreg <= 1/2.)                                                      //With probability 1/2, the individual will transmit type I meiotic product
    {
        vector<double> TransmChrom_LocPRDM9;                                //Determining the transmitted PRDM9 allele is the one of type I
        TransmChrom_LocPRDM9.push_back(IndivBuffer[0][0]);
        TransmChrom_LocPRDM9.push_back(IndivBuffer[0][1]);
        TransmittedChromosome.push_back(TransmChrom_LocPRDM9);              //PRDM9 is added to the transmitted meiotic product 2D table

        for(int i = 1; i < Nloc+1; ++i)
        {
            vector<double> TransmChrom_LocI;                                //Determining the transmitted target alleles are the one of type II
            TransmChrom_LocI.push_back(IndivBuffer[i][0]);
            TransmittedChromosome.push_back(TransmChrom_LocI);              //Each target is added to the transmitted meiotic product 2D table
        }
    }
    else                                                                    //With probability 1/2, the individual will transmit type II meiotic product
    {                                                                       //The equivalent algorithm to above is used to transmit type II
        vector<double> TransmChrom_LocPRDM9;
        TransmChrom_LocPRDM9.push_back(IndivBuffer[0][2]);
        TransmChrom_LocPRDM9.push_back(IndivBuffer[0][3]);
        TransmittedChromosome.push_back(TransmChrom_LocPRDM9);

        for(int i = 1; i < Nloc+1; ++i)
        {
            vector<double> TransmChrom_LocI;
            TransmChrom_LocI.push_back(IndivBuffer[i][1]);
            TransmittedChromosome.push_back(TransmChrom_LocI);
        }
    }

    double MutPRDM9 = unifDouble(mt);                                       //Drawing a double between 0 and 1

    if(MutPRDM9 < mu_A)                                                     //With probability mu_A, a mutation happens on the transmitted PRDM9 allele
    {
        TransmittedChromosome[0][1] = abs(TransmittedChromosome[0][1]-1);   //Allele A1 (0) mutates into A2 (1) and vice-versa
    }

    for(int k = 0; k < Nloc; ++k)                                           //At each target locus B
    {
        double Mut = unifDouble(mt);                                        //Drawing a double between 0 and 1

        if(Mut < mu_B)                                                      //With probability mu_B, a mutation happens on the transmitted target locus k+1
        {                                                                   //allele
            TransmittedChromosome[k+1][0] = abs(TransmittedChromosome[k+1][0]-1);       //Allele B1 (0) mutates into B2 (1) and vice-versa
        }
    }

    TransmChromAndFit.push_back(TransmittedChromosome);                     //The transmitted meiotic product 2D table is addes to the TransChromAndFit
                                                                            //3D table
    return TransmChromAndFit;
}

#endif // MEIOSIS_H_INCLUDED
