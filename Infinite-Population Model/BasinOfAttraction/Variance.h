#include <iostream>
#include <vector>
#include "Mean.h"
using namespace std;

long double Variance(vector<long double> X)
{
    int N = X.size();                   //N = number of elements in array X
    long double M = Mean(X);            //M = mean of X
    long double V(0.);

    for(int i = 0; i < N; ++i)
    {
        V += pow(X[i] - M,2);
    }

    V = V/N;                            //V = Sum(xi - M)^2/N : M = variance of X

    return V;
}
