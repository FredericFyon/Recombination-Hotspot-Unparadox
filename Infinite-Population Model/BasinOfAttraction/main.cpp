#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "Generation.h"
#include "ParameterConversion.h"
#include "Equilibre.h"
#include "EquilibreCorner.h"
#include "DistanceEq.h"
//#include "Mean.h"
#include "Variance.h"

//WITH MUT

using namespace std;

int main()
{                                                           //Parameters of the simulation:
    double b = 1.;                                          //b = probability that a binding results in a DSB
    double c = 1.;                                          //c = probability that a CO results in conversion
    double r_m = 0.5;                                         //r_m = recombination between PRDM9 and the target

    double MU[6] = {1.0e-12,1.0e-10,1.0e-8,1.0e-6,1.0e-5,2*1.0e-5};                 //MU = array of values of mu. mu = mutation rate

    /*vector<double> MU;
    for(int m = 0; m < 111; ++m)                            //Building MU with a large numbers of mu values
    {
        MU.push_back(pow(10,-15+m/10.));
    }*/


    double epsilon_int = 1.0e-2;                            //epsilon_int = small quantity with which compare the distance between a position and an internal equilibrium
    //double epsilon_corner = 1.0e-4;                       //epsilon_corner = small quantity with which compare the distance between a position and a corner equilibrium


/*                                                          //Creating a data output file for the average and variance of population mean fitness of every value of mu:
    std::ostringstream filemu;                              //An outside stringstream referred to as filemu is defined
    filemu << "/Users/fred/Dropbox/Fyon/RecombinationHotspost/Attraction Basin/AttractionBasinWithMutMAC/output_VarMuTest.dat";     //filemu is the path of the output file                              //This streams is made of a string defined here
    std::string varmu = filemu.str();                       //A string which is the streams content, defined above, is defined
    std::ofstream outMu(varmu);                             //An outside file, whose name is the string defined above, is defined
    outMu << std::setprecision(10);                         //The precision of the numerical values typed in the output file is set
*/


    for(int i = 0; i < 6; i++)                              //Results are calculated for the different mu values
    {
        double mu = MU[i];                                  //mu = mutation rate; it takes the ith value of the array MU
        cout << "mu = " << mu << " has begun" << endl;      //Display the current value of mu being evaluated

        //double f = 0.25;
        //double mu=1e-8;


/*
                                                            //Creating a data output file for ATTRACTION BASIN results of every mu for a given value of c
        std::ostringstream fileMut;                         //An outside stringstream referred to as fileMut is defined
        //fileMut << "/Users/fred/Dropbox/Fyon/RecombinationHotspost/Mutation Model/Basin_of_attraction/AttractionBasinWithMutMAC/Data_Basin_of_attraction/outputCycMut_c" << c << "_mu = " << mu << ".dat";  //fileMut is the path of the output file                                  //This streams is made of a string defined here
        fileMut << "/Users/Fred/Dropbox/Fyon/RecombinationHotspost/RecombDiv3_c" << c << "_mu = " << mu << ".dat";  //fileMut is the path of the output file                                  //This streams is made of a string defined here
        std::string varMut = fileMut.str();                 //A string which is the streams content, defined above, is defined
        std::ofstream outMut(varMut);                       //An outside file, whose name is the string defined above, is defined
        outMut << std::setprecision(4);                     //The precision of the numerical values typed in the output file is set


*/
                                                        //Creating a data output file for PERIOD results of every mu for a given value of c
        std::ostringstream fileCyc;                         //An outside stringstream referred to as fileCyc is defined
        fileCyc << "/Users/fred/Dropbox/Fyon/RecombinationHotspost/Mutation Model/Basin_of_attraction/AttractionBasinWithMutMAC/outputCycPeriod_mu = " << mu << ".dat";  //fileCyc is the path of the output file
        std::string varCyc = fileCyc.str();                 //A string which is the streams content, defined above, is defined
        std::ofstream outCyc(varCyc);                       //An outside file, whose name is the string defined above, is defined
        outCyc << std::setprecision(4);                   //The precision of the numerical values typed in the output file is set
        
                                                        //Creating a data output file for LIFE AND DEATH EXPECTANCIES results of every mu for a given value of c
        std::ostringstream fileLDE;                         //An outside stringstream referred to as fileCyc is defined
        fileLDE << "/Users/fred/Dropbox/Fyon/RecombinationHotspost/Mutation Model/Basin_of_attraction/AttractionBasinWithMutMAC/outputCycLifeAndDeath_mu = " << mu << ".dat";  //fileCyc is the path of the output file
        std::string varLDE = fileLDE.str();                 //A string which is the streams content, defined above, is defined
        std::ofstream outLDE(varLDE);                       //An outside file, whose name is the string defined above, is defined
        outLDE << std::setprecision(4);                   //The precision of the numerical values typed in the output file is set

        double f = 0.0001;                                    //f = fitness cost of no crossing-over; initiated at 0.001 because funny stuff happens at f = 0
        double d_f = 0.0001;                                  //df = incrementation of fitness cost: we will find the ATTRACTION BASIN / PERIOD for all values of f such as f = f + df, f < 1

        //vector<vector<double> > vecRec;
        //vector<double> f_array_Cycle;
        //vector<double> LogMu_array_Cycle;

        //double mu = 1e-10;
        //double d_mu = 0.001;

        while(f < c/2)                                    //Cycles only exist for f < c/2, there is no ATTRACTION BASIN nor PERIOD for f >= c/2
        //while(mu <= 1e-6)                                  //
        {
                                                             //Creating a data output file which compiles average and variance of mean fitnesses of population in cycling regime
           /* std::ostringstream fileCycW;                   //An outside stringstream referred to as fileCycW is defined
            fileCycW << "/Users/fred/Dropbox/Fyon/RecombinationHotspost/Attraction Basin/AttractionBasinWithMutMAC/outputCycW_mu = " << mu << "_f =" << f << ".dat";    //fileCycW is the path of the output file
            std::string varCycW = fileCycW.str();            //A string which is the streams content, defined above, is defined
            std::ofstream outCycW(varCycW);                  //An outside file, whose name is the string defined above, is defined
            outCycW << std::setprecision(6);                 //The precision of the numerical values typed in the output file is set

                                                             //Creating a data output file which compiles average and varianceof  mean fitnesses of population in internal equilibrium regime
            std::ostringstream fileEqW;                      //An outside stringstream referred to as fileEqW is defined
            fileEqW << "/Users/fred/Dropbox/Fyon/RecombinationHotspost/Attraction Basin/AttractionBasinWithMutMAC/outputEqW_mu = " << mu << "_f =" << f << ".dat";      //fileEqW is the path of the output file
            std::string varEqW = fileEqW.str();              //A string which is the streams content, defined above, is defined
            std::ofstream outEqW(varEqW);                    //An outside file, whose name is the string defined above, is defined
            outEqW << std::setprecision(6);*/                //The precision of the numerical values typed in the output file is set



            cout << "f = " << f << "..." << endl;            //Display the current value of f being evaluated
            //cout << "mu = " << mu << "..." << endl;            //Display the current value of mu being evaluated


            /*for(int j = 0; j < 41; ++j)                    //Calculation and iteration of mutation rates on locus A (PRDM9) and B (target)
            {
                double Ratio = pow(10,-2.+j/10.);            //Ratio = mua / mub. Results are calculated for different values of mua / mub
                double Eps = (Ratio-1.)/(1.+Ratio);
                double mua = MU[i]*(1+Eps);                  //mua = mutation rate on locus A
                double mub = MU[i]*(1-Eps);                  //mub = mutation rate on locus B
              */


                double mua = mu;                             //Calculation when mutation rates are the same between A and B:
                double mub = mu;                             //mua = mub = mu



                vector<double> newParameters = ParameterConversion(f,r_m,b,c);     //We define new parameters to simplify equations: alpha, beta, gamma, delta
                double alpha = newParameters[0];                                   //alpha is the first element of the array returned by ParameterConversion
                double beta = newParameters[1];                                    //beta is the second
                double gamma = newParameters[2];                                   //gamma is the third
                double delta = newParameters[3];                                   //delta is the fourth



                vector<double> X_eq = Equilibre(beta,gamma,delta,mu);              //We calculate the internal equilibrium, which depends on beta, gamma and delta
                double x1_eq = X_eq[0];                                            //A1B1 frequency at equilibrium
                double x2_eq = X_eq[1];                                            //A1B2 frequency at equilibrium
                double x3_eq = X_eq[2];                                            //A2B1 frequency at equilibrium
                double x4_eq = X_eq[3];



                //double d_x1 = 0.005;                                               //Incrementation of initial A1B1 frequencies: we want to find which frequencies lead to cycle and which to internal equilibrium
                double x1_init(0.005), x2_init(0.), x3_init(0.), x4_init(0.995);   //We define the initial frequencies first values. Initially, there is no A1B2 nor A2B1



                //while(x1_init < 1.)                                                //Initial frequency of A1B1 must be between 0 and 1
                //{
                    int t = 0;                                                     //t = current generation
                    int tmax = 2000000;                                            //We allow for a maximum of tmax generations


                    bool cycleFound = false;                                       //initially, we are not on the cycle, but the condition may become true at some point
                    bool internalEquilibrium = false;                              //initially, we are not in the internal equilibrium, but the condition may become true at some point
                    bool excessTime = false;                                       //excessTime will be true if we performed the tmax generation without reaching the cycle nor the internal equilibrium
                    bool xIncrease = false;                                        //xIncrease is true if the frequency of A1B1 increased compared to previous generation


                    long double x1(x1_init), x2(x2_init), x3(x3_init), x4(x4_init);//A1B1, A1B2, A2B1 and A2B2 are set at initial frequencies defined above
                    double Dist_Eq;                                                //Dist_Eq = distance to internal equilibrium


                    //vector<long double>WbarGen;                                  //WbarGen records mean population fitness every generation



                    vector<double> maxLocalGen, maxLocalX;                          //maxLocalGen records the generation time of every local maxima of frequency of A1B1
                                                                                    //maxLocalX records the frequency of A1B1 at every of those local maxima
            

                    while((cycleFound == false) && (internalEquilibrium == false) && (excessTime == false))   //Life cycle is performed until we reach the cycle, the internal equilibrium or the maximum generation time tmax
                    {
                        if(t == tmax)                                               //tmax is reached
                        {
                            excessTime = true;                                      //excessTime is true since we reached maximal generation time tmax without finding the cycle nor the internal equilibrium
                            cout << "No cycle nor internal equilibrium were found" << endl; //Display to the user that tmax was reached without finding the cycle nor the internal equilibrium
                        }



                        vector<long double> X = Generation(x1,x2,x3,x4,alpha,beta,gamma,delta,mua,mub);       //One generation is performed; X = array of haplotype frequencies after the generation
                        //WbarGen.push_back(X[4]);                                                            //The fifth element of X is the mean fitness of the population



                        if(abs(X[0] - x1_eq) < 0.1 && abs(X[1] - x2_eq) < 0.1 && abs(X[2] - x3_eq) < 0.1 && abs(X[3] - x4_eq) < 0.1) //Conditions: all 4 haplotype frequencies reasonably close to the ones of the internal equilibrium
                        {
                            Dist_Eq = DistanceEq(X[0],X[1],X[2],X[3],x1_eq,x2_eq,x3_eq,x4_eq);                //Dist_Eq = distance to internal equilibrium

                            if(Dist_Eq <= epsilon_int)                              //If Dist_Eq is small enough,
                            {
                                internalEquilibrium = true;                         //the internal equilibrium is considered to have been reached
                            }
                        }



                        if(X[0] > x1 + 1e-10)                                       //If frequency of A1B1 has increased compared to previous generation
                        {
                            xIncrease = true;                                       //it is recorded.
                        }
                        else if(X[0] < x1 - 1e-10)                                  //If frequency of A1B1 has decreased compared to previous generation
                        {
                            if(xIncrease == true)                                   //but xIncrease is true, which means that previously it had increased
                            {                                                       //a local maximum of frequency of A1B1 has been found.


                                if(maxLocalX.size() > 1 && (x1 - maxLocalX.back()) > -0.00004 && (x1 - maxLocalX.back()) < 0.00004)  //If the new local maximum is sufficiently close from the others, we are maybe in the limit cycle
                                {
                                    maxLocalX.push_back(x1);                        //so we add it the maximum to maxLocalX
                                    maxLocalGen.push_back(t+1);                     //and we add the current generation to maxLocalGen.



                                    if(maxLocalGen.size() > 19)                     //If there are already 20 local maxima or more recorded
                                    {                                               //we check if we have reached a stable cycle


                                        bool constantPeriod = true;                 //Let's say that the period between two local maxima is constant
                                        vector<int> diff;                           //diff = array of differences between two adjacent local maxima
                                        int k = maxLocalGen.size() - 10;            //k = number of local maxima - 10



                                        while(constantPeriod == true && k < maxLocalGen.size() - 1)     //We go through the last 10 local maxima
                                        {
                                            diff.push_back(maxLocalGen[k + 1] - maxLocalGen[k]);        //We add to diff the difference of generations between two current local maxima
                                            if(diff.size() > 1 && diff[k - maxLocalGen.size() + 10] - diff[k - maxLocalGen.size() + 9] > 10 && diff[k - maxLocalGen.size() + 10] - diff[k - maxLocalGen.size() + 9] < -10)
                                            {                                                           //If there are more than 10 generations of differences between two adjacent elements of 10
                                                constantPeriod = false;                                 //we consider that the cycle has not reached a stable period
                                            }
                                            k += 1;
                                        }



                                        if(k - maxLocalGen.size() + 1 == 0)                             //If we went through the last 10 local maxima without finding unstable period (constantPeriod = false)
                                        {
                                            cycleFound = true;                                          //then we consider that the limit cycle has been reached
                                            //outMut << f << " " << x1_init << endl;                    //We export the initial frequency of A1B1 in outMut: all initial frequencies recorded in outMut correspond to initial conditions leading to                                                             limit cycle: we thus get the ATTRACTION BASIN

                                          // if(x1_init <= 0.005 && j == 20)                            //If mua = mub and initial frequency is the smallest of the ones evaluated
                                            //{
                                                outCyc << f << " " << diff[k - maxLocalGen.size() + 9] << endl; //We export in outCyc the PERIOD of the cycle
                                           // }
                                        }
                                    }
                                }

                                else                                                //If the new local maximum is quite different from the others, we don't even bother checking if we are in the limit cycle
                                {
                                    maxLocalX.push_back(x1);                        //we record the local maxima, and carry on until finding the next one
                                }
                            }
                            xIncrease = false;                                      //We record that frequency of A1B1 has decreased
                        }


                        x1 = X[0];                                                  //We update x1, x2, x3 and x4 to prepare for next generation (unless cycle or equilibrium were found)
                        x2 = X[1];
                        x3 = X[2];
                        x4 = X[3];

                        t = t + 1;                                                  //One generation has been done, and we go to the next one (unless cycle or equilibrium were found)
                    }



                   /*                                                               //To record average and variance of population mean fitnesses:
                    if(internalEquilibrium == true)                                 //If the initial conditions ended in internal equilibrium
                    {
                        for(int p = 1; p < 10000; ++p)                              //10 000 supplementary generations are performed
                        {
                            vector<long double> X = Generation(x1, x2, x3, x4, alpha, beta, gamma, delta, mua, mub);
                            x1 = X[0];
                            x2 = X[1];
                            x3 = X[2];
                            x4 = X[3];
                            WbarGen.push_back(X[4]);                                //Mean population fitness is recorded for every generation
                        }


                        WbarGen.erase(WbarGen.begin(),WbarGen.end()-10000);         //WbarGen is made of the fitnesses during only those 10 000 last generations
                        long double WAmplitudeMean(Variance(WbarGen));              //Variance of population mean fitness is calculated
                        long double WbarMean(Mean(WbarGen));                        //Average of population mean fitness is calculated


                        outEqW << -2. + j/10.  << " " << WbarMean << " " << WAmplitudeMean << endl;     //In outEqW are recorded: the present mua / mub value, and the average and variance of population mean fitness
                        //outMu << log10(MU[i]) << " " << f << " " << WbarMean << " " << WAmplitudeMean << endl     //In outMu are recorded: the present mua = mub = mu value, and the average and variance of population mean fitness
                    }*/

                    if(cycleFound == true)                                          //If the initial conditions ended in limit cycle
                    {
                        
                        vector<double> X1,X2,X3,X4;                                     //Define four vectors that will record haplotype frequencies
                        
                        //f_array_Cycle.push_back(f);
                        //LogMu_array_Cycle.push_back(log10(mu));
                        vector<double> Rec;
                        for(int p = 1; p < 200000; ++p)                             //200 000 supplementary generations are performed
                        {
                            vector<long double> X = Generation(x1, x2, x3, x4, alpha, beta, gamma, delta, mua, mub);
                            x1 = X[0];
                            x2 = X[1];
                            x3 = X[2];
                            x4 = X[3];
                           // WbarGen.push_back(X[4]);                                //Mean population fitness is recorded for every generation
                            
                            
                            X1.push_back(x1);                                           //Record current haplotype frequencies
                            X2.push_back(x2);
                            X3.push_back(x3);
                            X4.push_back(x4);

                            Rec.push_back(r_m*b*(x1 + x4 - (x1*x4-x2*x3)));

                           /* if(p%20==0)
                            {
                                cout << p << endl;
                                outMut << x1 << " " << x2 << " " << x3 << " " << x4 << endl;
                            }*/

                        }
                        //vecRec.push_back(Rec);

                        /*WbarGen.erase(WbarGen.begin(),WbarGen.end()-200000);        //WbarGen is made of the fitnesses during only those 200 000 last generations
                        long double WAmplitudeMean(Variance(WbarGen));              //Variance of population mean fitness is calculated
                        long double WbarMean(Mean(WbarGen));                        //Average of population mean fitness is calculated
                        outCycW << -2. + j/10.  << " " << WbarMean << " " << WAmplitudeMean << endl;    //In outCycW are recorded: the present mua / mub value, and the average and variance of population mean fitness
                        //outMu << log10(MU[i]) << " " << f << " " << WbarMean << " " << WAmplitudeMean << endl;    //In outMu are recorded: the present mua = mub = mu value, and the average and variance of population mean fitness
                         */
                         //Calculation of the Life and Death Expectancies of Hotspots
                         
                         vector<double> HOT, COLD;                                                                                       //Vectors of generations in Hot and Cold states
                         
                         for(int tt = 0; tt < 200000; ++tt)                                               //Time spent in Hot or Cold is calculated in the time span that encompasses the last
                                                                                                                                         //10 local maxima (where we know the limit cycle is stable
                         {
                             
                             if(Rec[tt] >= 0.45)                                                                                          //If rate above 0.45, generation tt was in hot phase
                             {
                                 HOT.push_back(tt);
                             }
                             else if(Rec[tt]<= 0.05)                                                                                     //If rate below 0.05, generation tt was in cold phase
                             {
                                 COLD.push_back(tt);
                             }
                         }
                         
                         vector<double> StartHOT, EndHOT, StartCOLD, EndCOLD, LE, DE;                                                    //We record: start and end of hot phases, start and end of cold phases, life and death
                        double LifeExpectancy(0.),DeathExpectancy(0.);                                                               //expectancies
                         
                         int LLE(0),DDE(0);                                                                                              //Initialize duration of a given hot and cold phase at 0
                        
                         if(HOT.size() > 0)
                         {
                         for(int ttt = 0; ttt < HOT.size()-1; ++ttt)                                                                       //Running through the HOT vector
                         {
                             
                             if((HOT[ttt+1] == HOT[ttt] + 1) && StartHOT.size() > 0)                                                     //If two consecutive values of HOT are sequential generations, and if it's not the first
                                                                                                                                         //hot phase (which is truncated since we start at a local maximum of A1B1
                             {
                                 LLE = LLE+1;                                                                                            //We count the number of generations this happens. This will give the duration of one
                                                                                                                                         //Hot phase (one life expectancy)
                             }
                             else if(HOT[ttt] + 1 < HOT[ttt+1])                                                                          //If two consecutive values of HOT are NOT sequential generations, then it means there has
                                                                                                                                         //been a cold phase in between these generations
                             {
                                 if(LLE != 0)                                                                                            //If LLE was not 0, then we just run through a complete HOT phase, so we add LLE to the
                                                                                                                                         //vector of all HOT durations (vector of life expectancies)
                                 {
                                     LE.push_back(LLE);
                                     
                                     //cout << "HOT " << LLE << " " << StartHOT.back() << " " << HOT[ttt] << endl;
                                 }
                                 LLE = 0;                                                                                                //We reset LLE for the next HOT phase
                                 EndHOT.push_back(HOT[ttt]);                                                                             //We record that ttt has been the end of a HOT phase, and ttt+1 the beginning of a new one
                                 StartHOT.push_back(HOT[ttt+1]);
                             }
                         }
                                                                                          //Calculation of the average Life Expectancy
                             
                             for(int t4 = 0; t4 < LE.size();++t4)
                             {
                                 LifeExpectancy = LifeExpectancy + LE[t4];
                             }
                             
                             LifeExpectancy = LifeExpectancy/LE.size();
                             
                             
                         }
                        
                        if(COLD.size() > 0)
                        {
                         for(int ttt2 = 0; ttt2 < COLD.size()-1; ++ttt2)                                                                    //Same code applies for COLD phases and death expectancies
                         {
                             if((COLD[ttt2+1] == COLD[ttt2] + 1) && StartCOLD.size() > 0)
                             {
                                 DDE = DDE+1;
                             }
                             else if(COLD[ttt2] + 1 < COLD[ttt2+1])
                             {
                                 if(DDE != 0)
                                 {
                                     DE.push_back(DDE);
                                     
                                     //cout << "COLD " << DDE << " " << StartCOLD.back() << " " << COLD[ttt2] << endl;
                                 }
                                 DDE = 0;
                                 EndCOLD.push_back(COLD[ttt2]);
                                 StartCOLD.push_back(COLD[ttt2+1]);
                             }
                         }
                            
                            for(int t5 = 0; t5 < DE.size();++t5)
                            {
                                DeathExpectancy = DeathExpectancy + DE[t5];
                            }
                            
                            DeathExpectancy = DeathExpectancy/DE.size();
                        }
                         
                             outLDE << f << " " << LifeExpectancy << " " << DeathExpectancy << endl;
                  }

                    //x1_init = x1_init + d_x1;                                        //Once everything has been done and said, the whole thing is done for new initial frequencies of A1B1 and A2B2
                    //x2_init = 0.;
                    //x3_init = 0.;
                    //x4_init = 1. - x1_init;
                //}
            //}

            f = f + d_f;                                                             //Computations are done for various values of f
            //mu = pow(10,log10(mu) + d_mu);
        }
/*
        vector<double> Div;

        for(int j = 0; j < vecRec.size(); ++j)
        {
            vector<double> DiffRec;
            double Percent(0.);

            for(int k = 0; k < vecRec[j].size(); ++k)
            {
                DiffRec.push_back(abs(vecRec[j][k] - vecRec[(c/4-0.001-f_array_Cycle[0])/d_f + 1][k]));   //Comparison of recombination vectors between f = f0 + j and f = c/4
                //DiffRec.push_back(abs(vecRec[j][k] - vecRec[(-8-LogMu_array_Cycle[0])/d_mu][k]));
                if(DiffRec[k] < 0.03)
                {
                    Percent = Percent + 1;
                }
            }

            Percent = Percent/vecRec[j].size();

            Div.push_back(Percent);

            //outMut << (f_array_Cycle[0]+d_f*j) - c/4 + 0.001 << " " << Div[j] << endl;
            //outMut << (LogMu_array_Cycle[0]+d_mu*j) + 8 << " " << Div[j] << endl;
        }*/
    }                                                                                //Computations are done for various values of mu
}
