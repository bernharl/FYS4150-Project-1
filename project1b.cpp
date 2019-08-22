#include <iostream>
#include <fstream>
#include <cmath>
#include "armadillo"

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
    double a = -1; double b = 2; double c = -1;

    for (int i=0; i<=n-1;i++) 
    {
        for (int j=0; j<=n-1; j++)
        {
            if (i != j) {A[i][j] = 0;}
            
            else
            {
            A[i][j] = b;
            A[i][j+1] = c;
            A[i][j-1] = a;
            }
        cout << A[i][j];
        }
    cout<<endl;
    }
return 0;
}   