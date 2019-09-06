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
    int maxpow = 4; // Max power 10^maxpow
    double step = 0.5; // Step of power
    int arr_length = (int) (maxpow/step) - 1; // Length of arrays in loop
    mat maxerr = zeros <vec> (arr_length); // Array for saving maximum error
    mat time = zeros <vec> (arr_length); // Array for saving time taken
    mat n_arr = zeros <vec> (arr_length); // Array for saving power of 10
    mat numerical_solution;
    mat analytical_solution;
    for(float i=1; i <= maxpow; i+=step) // Looping over solutions with different power
    {
        int arrpos = (i - 1.0) / step; // Index of arrays to save values to
        int n = round(pow(10, i)); // n to use in an nxn matrix
        n_arr(arrpos) = n;
        double h = 1.0 / ((double) n + 1);

        double* x_arr = new double[n]; // x-values


        mat f_vec = zeros <vec> (n); // f values
        double h_sqr = h * h;
        for (int j = 0; j < n; j++) // Generating x and f-terms
        {
            x_arr[j] = (j + 1.0) * h;
            f_vec(j) = f(x_arr[j]) * h_sqr;
        }

        mat A = zeros <mat> (n, n); // Generating the tridiagonal matrix A using Armadillo
        A.diag(0).fill(2);
        A.diag(1).fill(-1);
        A.diag(-1).fill(-1);



        mat L, U;
        clock_t t_start = clock();
        lu(L, U, A); // LU decomoposition of A, saved to L and U

        mat y = solve(L, f_vec); // Solving intermediate step
        numerical_solution = solve(U, y); // Solving equation
        clock_t t_end = clock();
        time(arrpos) = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC; // Saving time spent calculating [ms]
        analytical_solution = zeros <vec> (n);
        for(int k=0; k<n; k++) // Calculating analytical solution for n points, used to get error
        {
            analytical_solution(k) = u(x_arr[k]);
        }
        mat tmp = log10(abs((numerical_solution - analytical_solution) / analytical_solution));
        maxerr(arrpos) =  tmp.max(); // Maximum log10 of absolute error for given ixi matrix solution
        delete [] x_arr;
    }
    // Saving vectors of interest
    n_arr.save("nLU.dat", csv_ascii);
    time.save("timeLU.dat", csv_ascii);
    maxerr.save("errorLU.dat", csv_ascii);
    numerical_solution.save("numLU.dat", csv_ascii);
    analytical_solution.save("anaLU.dat", csv_ascii);
    return 0;
}
