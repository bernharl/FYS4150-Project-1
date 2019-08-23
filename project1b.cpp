#include <iostream>
#include <fstream>
#include <cmath>
#include "armadillo"

using namespace std;

inline double f(double x)
{
    return 100 * exp(-10 * x);
}
inline double u(double x)
{
    return 1 - (1 - exp(-10)) * x - exp(-10*x);
}

int main(int argc, char *argv[]) 
{
    int n;
    string filename;

    filename = argv[1];
    n = atoi(argv[2]);
    double A[n][n];
    double a = -1; 
    double b = 2; 
    double c = -1;
    double h = 1.0/((double)n + 1);
    double h_sq = h*h;
    double v_arr[n];
    double b_arr[n];
    double x_arr[n];

    for (int i=1; i<=n; i++)
    {
        x_arr[i] = i*h;
        v_arr[i] = u(x_arr[i]);
        b_arr[i] = f(x_arr[i]) * h_sq;
        b_A[i] = b;
        //cout << x_arr[i] <<endl;
    }
    for (int i=1; i<=n-1; i++)
    {
        a_A[i] = a;
        c_A[i] = c;
    }
    

    return 0;
}
