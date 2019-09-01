#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
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

void data_to_file(double*, double*, double*, double*, string);
void Thomas_algorithm(string filename, int n, double &CPU_time_thomas, double &max_err_thomas, bool write = false)
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
    {   /* Filling up arrays for x-values, analytical function 
           values, r.h.s values and matrix diagonal */
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
    clock_t t_start = clock(); // Setting initial time for timing algorithm

    for (int i = 2; i < n; i++) 
    {   // Forward substitution
        b_A[i] = b_A[i] - a_A[i - 1] / b_A[i - 1] * c_A[i - 1];
        f_arr[i] = f_arr[i] - a_A[i - 1] / b_A[i - 1] * f_arr[i - 1];
    }

    v_arr[n - 1] = f_arr[n - 1]/b_A[n - 1]; // Solution at last point 

    for (int i = n - 2; i > 0; i--)
    {   // Backward substitution of remaining solution
        v_arr[i] = (f_arr[i] + v_arr[i + 1]) / b_A[i];
    }

    clock_t t_end = clock(); // End timing of algorithm
    CPU_time_thomas = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC; // Calculating CPU time [ms] of algorithm
    //cout << "CPU time: " << CPU_time << " ms " << endl; 


    for (int i = 1; i < n; i++)
    {   // Calculating relative error between analytical and numerical solution
        error[i] = fabs((u_arr[i] - v_arr[i]) / u_arr[i]);
    }

    if (write == true)
    {   /* Saving numerical values v, analytical calues u, x values 
           and log10 error in data file */
        data_to_file(v_arr, u_arr, x_arr, error, filename);
    }

    max_err_thomas = 0; 
    for (int i = 2; i < n; i++){
        // Calculating max error on x interval (0,1)
        if (error[i] > max_err_thomas){
            max_err_thomas = error[i];
        }
    }
    // Freeing allocated space in arrays
    delete [] a_A;
    delete [] b_A; 
    delete [] c_A;
    delete [] u_arr; 
    delete [] v_arr;
    delete [] f_arr; 
    delete [] x_arr;
    delete [] error;
}

void Specialized_algorithm(string filename, int n, double &CPU_time_special, double &max_err_special, bool write = false)
{
    //double a = -1.0;
    //double b = 2.0;

    double h = 1.0 / ((double) n);
    double h_sq = h * h;
    double* u_arr = new double[n + 1];
    double* v_arr = new double[n + 1];
    double* f_arr = new double[n + 1];
    double* x_arr = new double[n + 1];
    double* b_tilde = new double[n + 1];
    double* error = new double[n + 1];

    for (int i = 0; i <= n; i++)
    {
        x_arr[i] = i * h;
        u_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
    }
    v_arr[0] = v_arr[n] = 0;
    //algoritme
    b_tilde[0] = 2; b_tilde[n] = 2;
    for (int i = 1; i < n; i++)
    {
        b_tilde[i] = (i + 1.0) / ((double) i);
    }
    
    clock_t t_start = clock();

    for (int i = 2; i < n; i++)
    {
        f_arr[i] = f_arr[i] + f_arr[i - 1] / b_tilde[i - 1];
    }
    v_arr[n-1] = f_arr[n-1]/b_tilde[n-1];
    //cout << "v_arr[n - 1]" << v_arr[n - 1] << endl;
    for (int i = n - 2; i > 0; i--)
    {
        v_arr[i] = (f_arr[i] + v_arr[i + 1]) / b_tilde[i];
    }


    clock_t t_end = clock();
    CPU_time_special = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    //cout << "CPU time: " << CPU_time << " ms " << endl;
    
    for (int i = 1; i < n; i++)
    {
        //cout << u_arr[i] << "  " << f_arr[i] << endl;
        error[i] = fabs((u_arr[i] - v_arr[i]) / u_arr[i]);        
    }

    if (write == true){
        data_to_file(v_arr, u_arr, x_arr, error, filename);
    }

    max_err_special = 0;
    for (int i = 0; i < n; i++){
        if (error[i] > max_err_special){
            max_err_special = error[i];
        }
    }
    delete [] b_tilde; 
    delete [] u_arr; 
    delete [] v_arr;
    delete [] f_arr; 
    delete [] x_arr; 
    delete [] error;
}

void data_to_file(double *v_arr, double *u_arr, double *x_arr, double *error, string data_name)
{
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

void error_CPU_time_to_file(int exponent, string data_name, 
                            string error_name, string CPU_time_name)
{   
    ofstream outfile1;
    ofstream outfile2; 

    outfile1.open(error_name);
    outfile2.open(CPU_time_name);
    outfile1 << "n:" << "     log10(max_error)" << endl;
    outfile2 << "n:" << "     CPU time thomas [ms]:" << "        CPU time special [ms]:  " << endl;
    double CPU_time_special;
    double CPU_time_thomas;
    double max_err_special;
    double max_err_thomas;
    for (double i = 1; i <= exponent; i += 0.1)
    {   
        Thomas_algorithm(data_name, pow(10, i), CPU_time_thomas, max_err_thomas);
        Specialized_algorithm(data_name, pow(10, i), CPU_time_special, max_err_special);    
        outfile1 << pow(10, i) << " " << log10(max_err_special) << endl; 
        outfile2 << pow(10, i) << " " << CPU_time_special << " " << CPU_time_thomas << endl;
    }
    outfile1.close();
    outfile2.close();
}

void error_CPU_time_to_file( int, string, string, string);
int main()
{   
    error_CPU_time_to_file( exponent, data_name, error_name, CPU_time_name);
    //int n;
    //int exponent; 
    //string filename;

    //filename = argv[1];
    //n = atoi(argv[2]);
    //exponent = atoi(argv[3]);

    //double max_err, CPU_time = Thomas_algorithm(data_name, n, true);
    //double *r = Specialized_algorithm(data_name, n); 
    return 0;
}
