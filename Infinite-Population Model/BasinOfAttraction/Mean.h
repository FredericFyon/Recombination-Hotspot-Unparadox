#include <iostream>
#include <vector>
using namespace std;

long double Mean(vector<long double> X)
{
    int N = X.size();               //N = number of elements in array X
    long double M(0.);

    for(int i = 0; i < N; ++i)
    {
        M += X[i];
    }

    M = M/N;                        //M = Sum(xi)/N : M = mean of X

    return M;

}
