#include <iostream>
#include <vector>

vector<double> Equilibre(double beta, double gamma, double delta, double mu)
{
    double D;                               //D = linkage disequilibrium at internal equilibrium

    if(beta == gamma && mu == 0.)           //Special case for which D = 0
    {
        D = 0.;
    }
    else                                    //General case
    {
        double R = 4.*pow(delta + 2.*mu*(1.-delta),2) + pow(gamma,2)*pow(1.-2.*mu,2) - 2*beta*gamma*(1.-2.*mu) + pow(beta,2)*(1.-pow(2.*mu,2));
        D = (2*delta + 4*mu*(1. - delta) - sqrt(R))/(4.*(gamma - beta) - 8.*mu*(beta + gamma));             //Calculation of D following equation provided in Ubeda et al., 2021
    }

    double x1_eq = 0.25 + D;                //Calculation of haplotype frequencies at internal equilibrium following Ubeda et al., 2021
    double x2_eq = 0.25 - D;
    double x3_eq = 0.25 - D;
    double x4_eq = 0.25 + D;

    vector<double> X_eq;                    //Creation of the array composed of the 4 haplotype frequencies
    X_eq.push_back(x1_eq);
    X_eq.push_back(x2_eq);
    X_eq.push_back(x3_eq);
    X_eq.push_back(x4_eq);

    return X_eq;
}
