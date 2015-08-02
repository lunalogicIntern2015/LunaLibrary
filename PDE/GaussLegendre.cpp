#include <cmath>
#include "Quadrature.h"

//Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1]
 //and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula.

void gausslegendre(const double x1, const double x2, vector<double> &x, vector<double> &w) 
{
	const double EPS=1.0e-14;  
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	int n= (int)x.size(); 
	m=(n+1)/2;               
	xm=0.5*(x2+x1);            
	xl=0.5*(x2-x1);
	for (i=0;i<m;i++) {                        
		z=cos(3.141592654*(i+0.75)/(n+0.5));   
                                               
		do{                                    
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;j++)  
			{  
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);  
			}
			
			pp=n*(z*p1-p2)/(z*z-1.0);   
			z1=z;
			z=z1-p1/pp;                  
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;                   
		x[n-1-i]=xm+xl*z;               
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}

 /* 
//const double pi = 3.1415926535897932384626; 

//double func(const double x) 
//{
//  return 1/sqrt(2*pi)*exp(-x*x/2); // gaussian density~n(0,1) 
//}
//
//int main() 
//{ 
//	double x_min = -10; 
//	double x_max = 10; 
//    int N = 100; 
//    vector<double> x(N); 
//    vector<double> w(N); 
//    gausslegendre(x_min, x_max, x, w); 
//    double result = 0; 
//    for(int i=0; i<N; ++i) 
//	{
//	   result += func(x[i])*w[i]; 
//	}
// cout << "theoretical result = 1 vs numerical result " << result << endl; 
//}

 *///932384626; 

//double func(const double x) 
//{
//  return 1/sqrt(2*pi)*exp(-x*x/2); // gaussian density~n(0,1) 
//}
//
//int main() 
//{ 
//	double x_min = -10; 
//	double x_max = 10; 
//    int N = 100; 
//    vector<double> x(N); 
//    vector<double> w(N); 
//    gausslegendre(x_min, x_max, x, w); 
//    double result = 0; 
//    for(int i=0; i<N; ++i) 
//	{
//	   re