#pragma once   //! You should add this at the begin of each .h file 

//#include <cmath>
#include <vector>
//#include <PDE/NR/code/quadrature.h>
using namespace std;

//Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1]
 //and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula.


void gausslegendre(const double x1, const double x2, vector<double> &x, vector<double> &w); 
