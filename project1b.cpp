#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>


using namespace std;

inline double f(double x)
{
    return 100 * exp(-10 * x);
}
inline double u(double x)
{
    return 1 - (1 - exp(-10)) * x - exp(-10*x);
}

void Thomas_algorithm(string filename, int n)
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
    for (int i = 1; i <= n; i++) {
        //6 floating point operations
        //n - 1 iterations
        b_A[i] = b_A[i] - a_A[i-1] / b_A[i-1] * c_A[i-1];
        f_arr[i] = f_arr[i] - a_A[i-1] / b_A[i-1] * f_arr[i-1];

    }

    for (int i = n - 1; i >= 0; i--) {
        //3 floating point operations
        //n-1 iterations
        f_arr[i] = f_arr[i] - c_A[i] / b_A[i] * f_arr[i+1];
        //cout << f_arr[i]  << " i: " << i << endl;

    }
    //cout << "counter " << counter << endl;

    for (int i = 0; i < n; i++) {
        //1 floating point operation
        //n-1 iterations
        f_arr[i] /= b_A[i];

    }
    clock_t t_end = clock();
    double CPU_time = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    cout << "CPU time: " << CPU_time << " ms " << endl;
    /*
    ofstream outfile;
    outfile.open(filename);

    outfile << " v_i " << " u " << " x " << " Error " << endl;
    cout << u_arr[0] << endl;
    for (int i = 1; i < n; i++) {
        double error = abs((u_arr[i] - f_arr[i])/u_arr[i]);
        outfile << f_arr[i] << " " << u_arr[i] << " " << x_arr[i] << " "<< error << endl;
    }
    outfile.close();
    */
    int operations;
    operations = 6*(n - 1) + 3*(n - 1) + (n - 1);
    cout << operations << endl;

    delete [] a_A; delete [] b_A; delete [] c_A;
    delete [] u_arr; delete [] f_arr; delete [] x_arr;

}
void Specialized_algorithm(string filename, int n)
{
    double a = -1.0;
    double b = 2.0;
    double c = -1.0;

    double h = 1.0/((double) n + 1.0);
    double h_sq = h * h;
    double* u_arr = new double[n];
    double* f_arr = new double[n];
    double* x_arr = new double[n];

    for (int i = 0; i < n; i++)
    {
        x_arr[i] = i * h;
        u_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
        //cout << "Original " << f_arr[i] <<endl;
    }
    int counter = 0;


    clock_t t_start = clock();
    //algoritme
    double A = a - pow(b, 2)/a;
    double ba = b/a;
    double bA = b/A;

    for (int i=1; i<=n; i++)
    {
        f_arr[i] = f_arr[i] - ba * f_arr[i-1];
    }
    for (int i=n-1; i>=0; i--)
    {
        f_arr[i] = f_arr[i] - bA * f_arr[i+1];
    }
    for (int i = 0; i < n; i++) {
        //1 floating point operation
        //n-1 iterations
        //f_arr[i] /= b_A[i];
        f_arr[i] /= b;
    }

    /*
    ofstream outfile;
    outfile.open(filename);

    outfile << " v_i " << " u " << " x " << " Error " << endl;
    cout << u_arr[0] << endl;
    for (int i = 1; i < n; i++) {
        double error = abs((u_arr[i] - f_arr[i])/u_arr[i]);
        outfile << f_arr[i] << " " << u_arr[i] << " " << x_arr[i] << " "<< error << endl;
    }
    outfile.close();
    */
    clock_t t_end = clock();
    double CPU_time = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    cout << "CPU time: " << CPU_time << " ms " << endl;
    //delete [] a_A; delete [] b_A; delete [] c_A;
    delete [] u_arr; delete [] f_arr; delete [] x_arr;

}
int main(int argc, char *argv[])
{

    int n;
    string filename;

    filename = argv[1];
    n = atoi(argv[2]);

    Thomas_algorithm(filename, n);
    Specialized_algorithm(filename, n);

    return 0;
}
