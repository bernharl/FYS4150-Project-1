#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <armadillo>

using namespace arma;
using namespace std;

inline double f(double x)
{
    return 100.0 * exp(-10.0 * x);
}

int main(int argc, char *argv[])
{
    int n;
    string filename;

    filename = argv[1];
    n = atoi(argv[2]);    

    double h = 1.0 / ((double) n + 1.0);

    double* x_arr = new double[n];


    mat f_vec = zeros <mat> (n, 1);
    double h_sqr = h * h;
    for (int i = 0; i < n; i++)
    {
        x_arr[i] = (i + 1.0) * h;
        f_vec(i) = f(x_arr[i]) * h_sqr;
    }

    mat A = zeros <mat> (n, n);
    A.diag(0).fill(2);
    A.diag(1).fill(-1);
    A.diag(-1).fill(-1);
    
    

    mat L, U;
    lu(L, U, A);

    mat y = solve(L, f_vec);
    mat u = solve(U, y);


    cout << u << endl;
    u.save(filename, csv_ascii);
    return 0;
} 