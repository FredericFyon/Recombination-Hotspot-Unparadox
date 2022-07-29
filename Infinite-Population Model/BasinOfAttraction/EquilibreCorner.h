//
//  EquilibreInt.h
//  AttractionBasinWithMutMAC
//
//  Created by Frederic Fyon on 30/09/2020.
//

#ifndef EquilibreCorner_h
#define EquilibreCorner_h
#include <iostream>
#include <vector>

vector<double> EquilibreCorner(double beta, double gamma, double mu_a, double mu_b)
{
    double x1_eq = 1-(mu_b*beta+mu_a*(beta-gamma))*(1+beta)/(beta*(beta-gamma));
    double x2_eq = mu_b*(1+beta)/(beta-gamma);
    double x3_eq = mu_a*(1+beta)/beta;
    double x4_eq = 0.;

    vector<double> X_eq;
    X_eq.push_back(x1_eq);
    X_eq.push_back(x2_eq);
    X_eq.push_back(x3_eq);
    X_eq.push_back(x4_eq);

    return X_eq;
}


#endif /* EquilibreCorner_h */
