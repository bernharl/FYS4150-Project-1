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
    double A[n][n];
    double a = -1; 
    double b = 2; 
    double c = -1;
    double h = 1.0/((double)n + 1);
    double h_sq = h*h;
    double v_arr[n];
    double b_arr[n];
    double x_arr[n];

    for (int i=1; i<=n; i++)
    {
        x_arr[i] = i*h;
        v_arr[i] = u(x_arr[i]);
        b_arr[i] = f(x_arr[i]) * h_sq;
        //cout << x_arr[i] <<endl;
    }

    for (int i = 0; i <= n-1; i++){
        for (int j = 0; j<= n-1; j++){
            if (i == j) {A[i][j] = b;}
            else if (i == (j + 1)) {A[i][j] = c;}
            else if (i == j - 1) {A[i][j] = a;}
            else {A[i][j] = 0;}
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}
