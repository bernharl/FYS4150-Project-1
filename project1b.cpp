#include <iostream>
#include <fstream>
#include <cmath>

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
    double a = -1; 
    double b = 2; 
    double c = -1;
    double a_A[n-1];
    double b_A[n];
    double c_A[n-1];
    double h = 1.0/((double) n + 1);
    double h_sq = h*h;
    double v_arr[n];
    double f_arr[n];
    double x_arr[n];

    for (int i = 0; i < n; i++)
    {
        x_arr[i] = i*h;
        v_arr[i] = u(x_arr[i]);
        f_arr[i] = f(x_arr[i]) * h_sq;
        b_A[i] = b;
        cout << "Original " << f_arr[i] <<endl;
    }
    for (int i = 0; i < n - 1; i++)
    {
        a_A[i] = a;
        c_A[i] = c;
        
    }
    for (int i = 1; i <= n; i++){
        b_A[i] = b_A[i] - a_A[i-1] / b_A[i-1] * c_A[i-1];
        f_arr[i] = f_arr[i] - a_A[i-1] / b_A[i-1] * f_arr[i-1];
    }    

    for (int i = n - 1; i >= 0; i--){
        f_arr[i] = f_arr[i] - c_A[i] / b_A[i+1] * f_arr[i+1];
    }

    for (int i = 0; i < n; i++){
        f_arr[i] /= b_A[i];
    }

    for (int i = 0; i < n; i++){
        cout << f_arr[i]  << endl;
    }
    
    return 0;
}
