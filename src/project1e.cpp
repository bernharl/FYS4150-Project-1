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

inline double u(double x)
{   // Function returns analytical solution
    return 1.0 - (1.0 - exp(-10.0)) * x - exp(-10.0*x);
}


int main()
{
    int maxpow = 4;
    string filename_err = "errorLU.dat";
    mat maxerr = zeros <mat> (maxpow);

    for(int i=1; i <= maxpow; i++)
    {
        int n = pow(10, i);
        double h = 1.0 / ((double) n);

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
        mat numerical_solution = solve(U, y);

        mat analytical_solution = zeros <mat> (n);
        for(int ii=0; ii<n; ii++)
        {
            analytical_solution(ii) = u(x_arr[ii]);
        }
        mat tmp = log10(abs((numerical_solution - analytical_solution) / analytical_solution));
        maxerr(i-1) =  tmp.max();
    }

    maxerr.save(filename_err, csv_ascii);
    //numerical_solution.save("numtest.dat", csv_ascii);
    return 0;
}
