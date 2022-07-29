#include <iostream>
#include <cmath>

double DistanceEq(double x1, double x2, double x3, double x4, double x1_eq, double x2_eq, double x3_eq, double x4_eq)
{
    double Dist = sqrt(pow(x1 - x1_eq,2) + pow(x2 - x2_eq,2) + pow(x3 - x3_eq,2) + pow(x4 - x4_eq,2));
                    //Calculation of the Euclidean distance (4 dimensions) between the current frequencies (x1, x2, x3, x4)
                    //and the frequencies at internal equilibrium (x1_eq, x2_eq, x3_eq and x4_eq)


    return Dist;
}
