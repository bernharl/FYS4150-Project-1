#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace std;

inline double f(double x)
{
    return 100.0 * exp(-10.0 * x);
}
inline double u(double x)
{
    return 1.0 - (1.0 - exp(-10.0)) * x - exp(-10.0*x);
}

double Thomas_algorithm(string filename, int n, bool write = false)
{
    double a = -1.0;
    double b = 2.0;
    double c = -1.0;
    double h = 1.0/((double) n);
    double h_sq = h*h;
    double* a_A = new double[n];
    double* c_A = new double[n];
    double* b_A = new double[n + 1];
    double* u_arr = new double[n + 1];
    double* v_arr = new double[n + 1];
    double* f_arr = new double[n + 1];
    double* x_arr = new double[n + 1];
    double* error = new double[n + 1];

    for (int i = 0; i <= n; i++)
    {
        x_arr[i] = i * h;
        u_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
        b_A[i] = b;
    }
    v_arr[0] = v_arr[n] = 0;
    
    for (int i = 0; i < n; i++) {
        a_A[i] = a;
        c_A[i] = c;

    }
    clock_t t_start = clock();

    for (int i = 2; i < n; i++) 
    {
        b_A[i] = b_A[i] - a_A[i - 1] / b_A[i - 1] * c_A[i - 1];
        f_arr[i] = f_arr[i] - a_A[i - 1] / b_A[i - 1] * f_arr[i - 1];
    }

    v_arr[n - 1] = f_arr[n - 1]/b_A[n - 1];
    //cout << "v_arr[n - 1]" << v_arr[n - 1] << endl;
    for (int i = n - 2; i > 0; i--)
    {
        v_arr[i] = (f_arr[i] + v_arr[i + 1]) / b_A[i];
    }

    clock_t t_end = clock();
    double CPU_time = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    cout << "CPU time: " << CPU_time << " ms " << endl;


    for (int i = 1; i < n; i++)
    {   
        error[i] = fabs((u_arr[i] - v_arr[i]) / u_arr[i]);
    }

    if (write == true){
        ofstream outfile;
        outfile.open(filename);

        outfile << " v   " << "   u   " << "   x   " << "   log10(Error)   " << endl;
        //cout << u_arr[0] << endl;
        for (int i = 0; i <= n; i++) {
            outfile << v_arr[i] << " " << u_arr[i]
            << " " << x_arr[i] << " " << log10(error[i]) << endl;
        }

        outfile.close();
    }

    int operations;
    operations = 6 * (n - 1) + 3 * (n - 1) + (n - 1);
    //cout << operations << endl;
    double max_err = error[1];
    for (int i = 2; i < n; i++){
        if (error[i] > max_err){
            max_err = error[i];
        }
    }
    delete [] a_A; delete [] b_A; delete [] c_A;
    delete [] u_arr; delete [] f_arr; delete [] x_arr;
    delete [] error;
    return max_err;
}
double Specialized_algorithm(string filename, int n, bool write = false)
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
    double CPU_time = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    cout << "CPU time: " << CPU_time << " ms " << endl;
    
    for (int i = 1; i < n; i++)
    {
        //cout << u_arr[i] << "  " << f_arr[i] << endl;
        error[i] = fabs((u_arr[i] - v_arr[i]) / u_arr[i]);        
    }

    if (write == true){
        ofstream outfile;
        outfile.open(filename);

        outfile << " v   " << "   u   " << "   x   " << "   log10(Error)   " << endl;
        //cout << u_arr[0] << endl;
        for (int i = 0; i <= n; i++) {
            outfile << v_arr[i] << " " << u_arr[i]
            << " " << x_arr[i] << " " << log10(error[i]) << endl;
        }

        outfile.close();
    }
    int max_i   = 0;
    double max_err = 0;
    for (int i = 0; i < n; i++){
        if (error[i] > max_err){
            max_err = error[i];
            max_i = i;
        }
    }
    cout << max_i << endl;
    delete [] b_tilde; 
    delete [] u_arr; 
    delete [] v_arr;
    delete [] f_arr; 
    delete [] x_arr; 
    delete [] error;
    
    return log10(max_err);
}
int main(int argc, char *argv[])
{

    int n;
    int exponent; 
    string filename;

    filename = argv[1];
    n = atoi(argv[2]);
    exponent = atoi(argv[3]);

    //double max_err = Thomas_algorithm(filename, n, true);
    //double max_err = Specialized_algorithm(filename, n, true);
    //cout << max_err << endl;
    
    ofstream outfile;
    outfile.open("error.dat");
    outfile << "n" << "     log10(max_error)" << endl;
    for (double i = 1; i <= exponent; i += 0.1)
    {
        outfile << pow(10, i) << " " << Specialized_algorithm(filename, pow(10, i)) << endl;
    }
    
    outfile.close();
    
    return 0;
}
