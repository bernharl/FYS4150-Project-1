#include <iostream>
#include <string>
using namespace std;

// Single grid size 
int n = 1e2;        
// Upper level of exponents to run over in max(error) vs n
int exponent = 7;   
// Name of data file used for plotting u[i] and v[i] vs x[i]
string data_name = "data.dat";  
// Name of thomas file to save max log10(error) and CPU time vs n
string thomas_name = "thomas.dat";  
// Name of specialised thomas file to save max log10(error) and CPU time vs n
string special_name = "special.dat";    
// Name of data file used for specialised Thomas algorithm 
string data_name_special = "data_special.dat";
// Name of data file used for Thomas algorithm 
string data_name_thomas = "data_thomas.dat";
