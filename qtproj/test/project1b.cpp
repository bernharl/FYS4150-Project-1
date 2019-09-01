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
    double* a_A = new double[n-1];

    double* c_A = new double[n-1];

    double* b_A = new double[n];

    double h = 1.0/((double) n + 1.0);
    double h_sq = h*h;
    double* u_arr = new double[n];
    double* f_arr = new double[n];
    double* x_arr = new double[n];
    double* error = new double[n];

    for (int i = 0; i < n; i++)
    {
        x_arr[i] = i*h;
        u_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
        b_A[i] = b;
        //cout << "Original " << f_arr[i] <<endl;
    }
    int counter = 0;
    for (int i = 0; i < n - 1; i++) {
        a_A[i] = a;
        c_A[i] = c;

    }
    clock_t t_start = clock();
    for (int i = 1; i < n; i++) {
        //6 floating point operations
        //n - 1 iterations
        b_A[i] = b_A[i] - a_A[i-1] / b_A[i-1] * c_A[i-1];
        f_arr[i] = f_arr[i] - a_A[i-1] / b_A[i-1] * f_arr[i-1];
        //cout << f_arr[i] << endl;

    }

    for (int i = n-2; i > 0; i--) {
        //3 floating point operations
        //n-1 iterations

        f_arr[i] = f_arr[i] - c_A[i] / b_A[i] * f_arr[i+1];
        //cout << f_arr[i] << endl;


    }
    //cout << "counter " << counter << endl;

    for (int i = 0; i < n; i++) {
        //1 floating point operation
        //n-1 iterations
        f_arr[i] /= b_A[i];
        //cout << f_arr[i] << endl;

    }
    clock_t t_end = clock();
    double CPU_time = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    cout << "CPU time: " << CPU_time << " ms " << endl;


    for (int i = 1; i < n; i++){
        error[i] = abs((u_arr[i] - f_arr[i]) / u_arr[i]);
    }
    if (write == true){
        ofstream outfile;
        outfile.open(filename);

        outfile << " v_i   " << "   u   " << "   x   " << "   Error   " << endl;
        //cout << u_arr[0] << endl;
        for (int i = 1; i < n; i++) {
            outfile << f_arr[i] << " " << u_arr[i]
            << " " << x_arr[i] << " " << error[i] << endl;
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
    double b = 2.0;

    double h = 1.0 / ((double) n + 1.0);
    double h_sq = h * h;
    double* u_arr = new double[n];
    double* f_arr = new double[n];
    double* x_arr = new double[n];
    double* b_tilde = new double[n];
    double* error = new double[n];

    for (int i = 0; i < n; i++)
    {
        x_arr[i] = i * h;
        u_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
        //cout << "Original " << f_arr[i] <<endl;
        b_tilde[i] = 2.0 - i / (i + 1.0);
        //cout << b_tilde[i] << "  " << i << endl;

    }

    int counter = 0;

    //algoritme
    //double b_tilde = b - pow(a, 2) / b;
    b_tilde[0] = 2;
    for (int i=1; i<=n; i++)
    {
      //b_tilde[i] = 2.0 - (i - 1.0) / i;
      b_tilde[i] = 1.0 + 1.0 / i;
      //b_tilde[i] = b_tilde[i] - 1.0 / b_tilde[i-1];
      //cout << b_tilde[i] << endl;

    }
    //cout << A << endl;
    //double ab = -1.0 / b;
    //cout << ba << endl;
    //double ab_tilde = -1.0 / b_tilde;
    //cout << bA << endl;
    clock_t t_start = clock();
    for (int i=1; i<n; i++)
    {
        //f_arr[i] = f_arr[i] - c_A[i] / b_A[i] * f_arr[i+1];
        f_arr[i] = f_arr[i] + f_arr[i-1] / b_tilde[i];
        //cout << f_arr[i] << endl;
    }

    for (int i=n-2; i>0; i--)
    {

        f_arr[i] = f_arr[i] + f_arr[i+1] / b_tilde[i+1];
        //cout << f_arr[i] << endl;

    }

    for (int i = 0; i < n; i++) {
        //1 floating point operation
        //n-1 iterations
        //f_arr[i] /= b_A[i];
        f_arr[i] = f_arr[i] / b_tilde[i+1];
        //cout << f_arr[i] << endl;
    }
    clock_t t_end = clock();
    double CPU_time = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    cout << "CPU time: " << CPU_time << " ms " << endl;

    for (int i = 1; i < n; i++){
        error[i] = abs((u_arr[i] - f_arr[i]) / u_arr[i]);
    }
    if (write == true){
        ofstream outfile;
        outfile.open(filename);

        outfile << " v_i   " << "   u   " << "   x   " << "   Error   " << endl;
        //cout << u_arr[0] << endl;
        for (int i = 1; i < n; i++) {
            outfile << f_arr[i] << " " << u_arr[i]
            << " " << x_arr[i] << " " << error[i] << endl;
        }

        outfile.close();
    }
    double max_err = error[1];
    for (int i = 2; i < n; i++){
        if (error[i] > max_err){
            max_err = error[i];
        }
    }

    cout << "jeff" << endl;
    for (int i = 0; i < n-1; i++){
        cout << b_tilde[i] << endl;
    }
    cout << "jeff" << endl;
    cout << "jeff" << endl;
    delete [] u_arr; 
    cout << "jeff" << endl;
    delete [] f_arr; 
    cout << "jeff" << endl;

    delete [] x_arr; 
    cout << "jeff" << endl;
    delete [] error;
    cout << "jeff" << endl;

    return max_err;
}
int main(int argc, char *argv[])
{

    int n;
    string filename;

    filename = argv[1];
    n = atoi(argv[2]);

    //double max_err = Thomas_algorithm(filename, n);
    double max_err = Specialized_algorithm(filename, 201);
    cout << max_err << endl;
    //ofstream outfile;
    //outfile.open("error.dat");
    //outfile << "n" << "max_error" << endl;
    //outfile << 201 << " " << Specialized_algorithm(filename, 201) << endl;
    /*
    for (int i = 1; i < 1000000; i += 100){
        cout << i << endl;
    }
    */
    //outfile.close();
    return 0;
}