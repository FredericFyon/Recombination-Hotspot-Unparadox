#include <iostream>
#include <vector>
using namespace std;

vector<long double> Generation(long double x1, long double x2, long double x3, long double x4, double alpha, double beta, double gamma, double delta, double mu_a, double mu_b)
{
    long double y_sel_1, y_sel_2, y_sel_3, y_sel_4, y_mut_1, y_mut_2,y_mut_3,y_mut_4;       //y_sel = frequencies after selection, y_mut = frequencies after mutation
    long double D = x1*x4 - x2*x3;                                                          //D = linkage disequilibrium at current generation

    y_sel_1 = ((alpha + beta*x1 - gamma*x2)*x1 - delta*D);                                  //Calculation of frequencies after selection, following equation provided in Ubeda et al. 2021
    y_sel_2 = ((alpha - beta*x2 + gamma*x1)*x2 + delta*D);
    y_sel_3 = ((alpha - beta*x3 + gamma*x4)*x3 + delta*D);
    y_sel_4 = ((alpha + beta*x4 - gamma*x3)*x4 - delta*D);

    long double Wbar = y_sel_1 + y_sel_2 + y_sel_3 + y_sel_4;                          //Wbar = population mean fitness

    y_sel_1 = y_sel_1/Wbar;                                                            //Normalization of frequencies after selection
    y_sel_2 = y_sel_2/Wbar;
    y_sel_3 = y_sel_3/Wbar;
    y_sel_4 = y_sel_4/Wbar;

    y_mut_1 = (1 - mu_a - mu_b + mu_a*mu_b)*y_sel_1 + mu_b*(1 - mu_a)*y_sel_2 + mu_a*(1 - mu_b)*y_sel_3 + mu_a*mu_b*y_sel_4;          //Calculation of frequencies after mutation
    y_mut_2 = (1 - mu_a - mu_b + mu_a*mu_b)*y_sel_2 + mu_b*(1 - mu_a)*y_sel_1 + mu_a*(1 - mu_b)*y_sel_4 + mu_a*mu_b*y_sel_3;
    y_mut_3 = (1 - mu_a - mu_b + mu_a*mu_b)*y_sel_3 + mu_b*(1 - mu_a)*y_sel_4 + mu_a*(1 - mu_b)*y_sel_1 + mu_a*mu_b*y_sel_2;
    y_mut_4 = (1 - mu_a - mu_b + mu_a*mu_b)*y_sel_4 + mu_b*(1 - mu_a)*y_sel_3 + mu_a*(1 - mu_b)*y_sel_2 + mu_a*mu_b*y_sel_1;


    vector<long double> newGeneration;                                                      //Creation of the array containing the new haplotype frequencies after one generation
    newGeneration.push_back(y_mut_1);
    newGeneration.push_back(y_mut_2);
    newGeneration.push_back(y_mut_3);
    newGeneration.push_back(y_mut_4);
    newGeneration.push_back(Wbar);                                                          //Fifth element of the array is the population mean fitness

    return newGeneration;
}
