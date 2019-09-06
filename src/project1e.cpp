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
    double step = 0.5;
    int arr_length = (int) (maxpow/step) - 1;
    string filename_err = "errorLU.dat";
    mat maxerr = zeros <mat> (arr_length);
    mat numerical_solution;
    mat analytical_solution;
    for(float i=1; i <= maxpow; i+=step)
    {
        int arrpos = (i - 1.0) / step;
        cout << arrpos << endl;
        int n = round(pow(10, i));
        double h = 1.0 / ((double) n + 1);

        double* x_arr = new double[n];


        mat f_vec = zeros <vec> (n);
        double h_sqr = h * h;
        for (int j = 0; j < n; j++)
        {
            x_arr[j] = (j + 1.0) * h;
            f_vec(j) = f(x_arr[j]) * h_sqr;
        }

        mat A = zeros <mat> (n, n);
        A.diag(0).fill(2);
        A.diag(1).fill(-1);
        A.diag(-1).fill(-1);



        mat L, U;
        lu(L, U, A);

        mat y = solve(L, f_vec);
        numerical_solution = solve(U, y);

        analytical_solution = zeros <vec> (n);
        for(int k=0; k<n; k++)
        {
            analytical_solution(k) = u(x_arr[k]);
        }
        mat tmp = log10(abs((numerical_solution - analytical_solution) / analytical_solution));
        maxerr(arrpos) =  tmp.max();
        cout << tmp.index_max() << " " << tmp.max() << endl;
        delete [] x_arr;
    }
    maxerr.print();
    maxerr.save(filename_err, csv_ascii);
    numerical_solution.save("numLU.dat", csv_ascii);
    analytical_solution.save("anaLU.dat", csv_ascii);
    // cout << maxerr << endl;
    return 0;
}
