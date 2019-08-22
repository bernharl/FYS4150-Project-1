#include <iostream>
#include <fstream>
#include <cmath>
#include "armadillo"

using namespace std;

inline double f(double x)
{
    return 100 * exp(-10 * x)
};
inline double u(double x)
{
    return 1 - (1 - exp(-10)) * x - exp(-10*x)
};

int main(int argc, char *argv[]) 
{
    int dimension;
    string filename;
    
    
}