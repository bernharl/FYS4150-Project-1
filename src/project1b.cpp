#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include "keyword_args.h"

using namespace std;

inline double f(double x)
{   /* Function returns right side of matrix equation
       A*v = f */
    return 100.0 * exp(-10.0 * x);
}

inline double u(double x)
{   // Function returns analytical solution
    return 1.0 - (1.0 - exp(-10.0)) * x - exp(-10.0*x);
}

void data_to_file(double*, int, double*, double*, double*, string);
void Thomas_algorithm(string filename, int n, double &CPU_time_thomas,
                        double &max_err_thomas, bool write = false)
{   /* Function solves one-dimensional Poisson equation
       using a Thomas algorithm. Function takes filename
       of file to save data to, matrix/array dimensions n
       and a boolean argument write to determain if data
       file is produced. */
    double a = -1.0;    // Upper diagonal elements
    double b = 2.0;     // Diagonal elements
    double c = -1.0;    // Lower diagonal elements
    double h = 1.0/((double) n);    // Step size
    double h_sq = h*h;
    double* a_A = new double[n];    // Defining upper matrix diagonal as array
    double* c_A = new double[n];    // Defining matrix diagonal as array
    double* b_A = new double[n + 1];    // Defining lower matrix diagonal as array
    double* u_arr = new double[n + 1];  // Array for analytical solution
    double* v_arr = new double[n + 1];  // Array for numerical soluton
    double* f_arr = new double[n + 1];  // Array for r.h.s of the matrix equation
    double* x_arr = new double[n + 1];  // Array for x-axis
    double* error = new double[n + 1];  // Array for relative errors as function of x

    for (int i = 0; i <= n; i++)
    {   /* Filling up arrays for x-values, analytical solution,
           r.h.s values and matrix diagonal */
        x_arr[i] = i * h;
        u_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
        b_A[i] = b;
    }
    v_arr[0] = v_arr[n] = 0; // Setting boundary condition for numerical solution

    for (int i = 0; i < n; i++)
    {   // Filling up upper and lower matrix diagonals
        a_A[i] = a;
        c_A[i] = c;

    }
    clock_t t_start = clock(); // Initializing timer

    for (int i = 2; i < n; i++)
    {   // Forward substitution
        b_A[i] = b_A[i] - a_A[i - 1] / b_A[i - 1] * c_A[i - 1];
        f_arr[i] = f_arr[i] - a_A[i - 1] / b_A[i - 1] * f_arr[i - 1];
    }

    v_arr[n - 1] = f_arr[n - 1]/b_A[n - 1]; // Solution at last grid point

    for (int i = n - 2; i > 0; i--)
    {   // Backward substitution to find remaining solutions
        v_arr[i] = (f_arr[i] - c_A[i] * v_arr[i + 1]) / b_A[i];
    }

    clock_t t_end = clock(); // End timer
    CPU_time_thomas = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC; // Calculating CPU time [ms]


    for (int i = 1; i < n; i++)
    {   // Calculating relative error between analytical and numerical solution
        error[i] = fabs((u_arr[i] - v_arr[i]) / u_arr[i]);
    }

    if (write == true)
    {   /* Saving numerical values v, analytical calues u, x values
           and log10 error in data file */
        data_to_file(v_arr, n, u_arr, x_arr, error, filename);
    }

    max_err_thomas = 0; // Initializing max relative error
    for (int i = 2; i < n; i++)
    {   // Calculating max error on x interval (0,1)
        if (error[i] > max_err_thomas){
            max_err_thomas = error[i];
        }
    }
    // Freeing allocated memory
    delete [] a_A;
    delete [] b_A;
    delete [] c_A;
    delete [] u_arr;
    delete [] v_arr;
    delete [] f_arr;
    delete [] x_arr;
    delete [] error;
}

void Specialized_algorithm(string filename, int n, double& CPU_time_special,
                            double& max_err_special, bool write = false)
{   /* Function solves one-dimensional Poisson equation
       using a specialized Thomas algorithm. Function takes filename
       of file to save data to, matrix/array dimensions n,
       CPU time variable and max error variables to compute
       and a boolean argument write to determain if data
       file is produced. */

    double h = 1.0 / ((double) n);  // Step size
    double h_sq = h * h;
    double* u_arr = new double[n + 1];  // Arraay for analytical solution
    double* v_arr = new double[n + 1];  // Array for numerical solution
    double* f_arr = new double[n + 1];  // Array for r.h.s of matrix equation
    double* x_arr = new double[n + 1];  // Array for x-axis
    double* b_tilde = new double[n + 1];    // Array for updated diagonal elements
    double* error = new double[n + 1];      // Arrray for relative error

    for (int i = 0; i <= n; i++)
    {   /* Filling up arrays for x-values, analytical solution,
         and r.h.s of matrix equation */
        x_arr[i] = i * h;
        u_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
    }
    v_arr[0] = v_arr[n] = 0; // Setting boundary conditions for numerical solution

    b_tilde[0] = 2; b_tilde[n] = 2; // Setting endpoints of diagonal vector
    for (int i = 1; i < n; i++)
    {   //Updating diagonal elements due to Gauss elimination
        b_tilde[i] = (i + 1.0) / ((double) i);
    }

    clock_t t_start = clock();  // Initializing timer

    for (int i = 2; i < n; i++)
    {   // Forward substitution
        f_arr[i] = f_arr[i] + f_arr[i - 1] / b_tilde[i - 1];
    }

    v_arr[n-1] = f_arr[n-1]/b_tilde[n-1]; // Solution at last grid point
    for (int i = n - 2; i > 0; i--)
    {   // Backward substitution to find remaining solutions
        v_arr[i] = (f_arr[i] + v_arr[i + 1]) / b_tilde[i];
    }


    clock_t t_end = clock(); // Ending timer
    CPU_time_special = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC; // Calculating CPU time [ms]


    for (int i = 1; i < n; i++)
    {   // Calculating relative error between analytical and numerical solution
        error[i] = fabs((u_arr[i] - v_arr[i]) / u_arr[i]);
    }

    if (write == true)
    {   /* Saving numerical values v, analytical calues u, x values
           and log10 error in data file */
        data_to_file(v_arr, n, u_arr, x_arr, error, filename);
    }

    max_err_special = 0;    // Initializing max relative error
    for (int i = 0; i < n; i++)
    {   // Calculating max error on x interval (0,1)
        if (error[i] > max_err_special){
            max_err_special = error[i];
        }
    }
    // Freeing allocated memory
    delete [] b_tilde;
    delete [] u_arr;
    delete [] v_arr;
    delete [] f_arr;
    delete [] x_arr;
    delete [] error;
}

void data_to_file(double *v_arr, int n, double *u_arr, double *x_arr, double *error, string data_name)
{   /* Function prints numerical solution, analytical solution,
       x values and log10 of the relative error to
       data file for single grid size n.*/

    ofstream outfile;
    outfile.open(data_name);
    outfile << " v   " << "   u   " << "   x   " << "   log10(Error)   " << endl;
    for (int i = 0; i <= n; i++)
    {
        outfile << v_arr[i] << " " << u_arr[i]
        << " " << x_arr[i] << " " << log10(error[i]) << endl;
    }
    outfile.close();
}

void thomas_n_to_file(int exponent, string data_name, string thomas_name)
{   /* Function prints log10 of relative error and CPU time
       for Thomas algorithm for different grid sizes n */
    ofstream outfile;
    outfile.open(thomas_name);
    outfile << "n:" << setw(20) <<  "log10(Max error):" << setw(20) << "CPU time [ms]:" << endl;
    double CPU_time_thomas;
    double max_err_thomas;
    for (double i = 1; i <= exponent; i += 0.1)
    {
        Thomas_algorithm(data_name, round(pow(10, i)), CPU_time_thomas, max_err_thomas);
        outfile << round(pow(10, i)) << setprecision(10) << setw(20) << log10(max_err_thomas)
         << setprecision(10) << setw(20) << CPU_time_thomas << setprecision(10) << setw(20) << endl;
    }
    outfile.close();
}

void special_n_to_file(int exponent, string data_name, string special_name)
{   /* Function prints log10 of relative error and CPU time
       for specialised thomas algorithm for different grid sizes n */
    ofstream outfile;
    outfile.open(special_name);
    outfile << "n:" << setw(20) <<  "Max error:" << setw(20) << "CPU time [ms]:" << endl;
    double max_err_special;
    double CPU_time_special;
    for (double i = 1; i <= exponent; i += 0.1)
    {
        Specialized_algorithm(data_name, round(pow(10, i)), CPU_time_special, max_err_special);
        outfile << round(pow(10, i)) << setprecision(10) << setw(20) << log10(max_err_special)
        << setprecision(10) << setw(20) << CPU_time_special << setprecision(10) << setw(20) << endl;
    }
    outfile.close();
}

int main()
{
    thomas_n_to_file(exponent, data_name, thomas_name);
    special_n_to_file(exponent, data_name, special_name);

    double dummyx;
    double dummyxx;
    Specialized_algorithm(data_name_special, 10000, dummyx, dummyxx, true);
    Thomas_algorithm(data_name_thomas, 10000, dummyx, dummyxx, true);

    //int n;
    //int exponent;
    //string filename;

    //filename = argv[1];
    //n = atoi(argv[2]);
    //exponent = atoi(argv[3]);
    //double max_err_thomas;
    //double CPU_time_thomas;
    //Thomas_algorithm(data_name, n, true, CPU_time_thomas, max_err_thomas);
    //double *r = Specialized_algorithm(data_name, n);
    return 0;
}
