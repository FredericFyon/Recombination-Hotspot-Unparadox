#include <iostream>
#include <vector>
using namespace std;

vector<double> ParameterConversion(double f, double r_m, double b, double c)
{
    double alpha, beta, gamma, delta;                                   //alpha, beta, gamma, delta are the new parameters

    alpha = 0.25*b + 0.5*(1.-0.5*b)*(1.-f);                             //Calculation of alpha following Ubeda et al., 2021
    beta = 0.25*b*f;                                                    //Calculation of beta following Ubeda et al., 2021
    gamma = 0.125*b*c;                                                  //Calcultion of gamma following Ubeda et al., 2021
    delta = 0.125*b*(0.5*c + (1.-c)*r_m)+0.25*r_m*(1.-0.5*b)*(1.-f);    //Calculation of delta following Ubeda et al., 2021


    vector<double> newParameters;                                       //Creation of the array containing the new parameters alpha, beta, gamma, delta
    newParameters.push_back(alpha);
    newParameters.push_back(beta);
    newParameters.push_back(gamma);
    newParameters.push_back(delta);

    return newParameters;
}
