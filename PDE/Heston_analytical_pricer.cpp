#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "Interpolator.h"
#include "Heston_analytical_pricer.h"
#include "Params.h"
#include "Quadrature.h"

using namespace std;

//! second forumula for the characteristic Heston: "The little Hestion trap" : https://perswww.kuleuven.be/~u0009713/HestonTrap.pdf

//! constant constant
const double pi_inverse = 1.0/3.1415926535897932384626;
const complex<double> complex_i = complex<double>(0,1);
const complex<double> complex_unity = complex<double>(1,0);

bool OPTIMIZATION_if_use_CF_second_formula = true;  //! make the v==0 case Heston's analytical price don't behave like sin! 

Heston_analytical_pricer :: Heston_analytical_pricer()
{
	      //! Hestion's Params
		  kappa_   = Heston_params_kappa;
		  theta_   = Heston_params_theta;
		  sigma_   = Heston_params_sigma;
		  rho_     = Heston_params_rho;
		  lambda_  = 0.0;    // corresponding to different risk-neutral probability
		  r_       = Heston_params_r;
		  a_       = Heston_params_kappa*Heston_params_theta;

		  //! Option's params
		  K_ = Option_params_K;
		  T_ = Option_params_T;

		  //! Numerical Integral's Params
		  integral_interval_begin_ = 0.00000001;
		  integral_interval_end_   = 100.0;
}

Heston_analytical_pricer :: Heston_analytical_pricer(double K)
{
	      //! Hestion's Params
		  kappa_   = Heston_params_kappa;
		  theta_   = Heston_params_theta;
		  sigma_   = Heston_params_sigma;
		  rho_     = Heston_params_rho;
		  lambda_  = 0.0;    // corresponding to different risk-neutral probability
		  r_       = Heston_params_r;
		  a_       = Heston_params_kappa*Heston_params_theta;

		  //! Option's params
		  K_ = K;
		  T_ = Option_params_T;

		  //! Numerical Integral's Params
		  integral_interval_begin_ = 0.00000001;
		  integral_interval_end_   = 100.0;
}

////! don't care the efficiency, firstly accuracy.
////! OK
//complex<double> Heston_analytical_pricer :: f(double S, double v,
//											  double b, double u,
//											  double phi)
//{
//
//	double tau = T_;
//	double x = log(S);
//
//	complex<double> d = sqrt( 
//					 (rho_*sigma_*phi*complex_i-b)*(rho_*sigma_*phi*complex_i-b) - sigma_*sigma_*(2*u*phi*complex_i - phi*phi)
//				    );
//
//	complex<double> g = (b-rho_*sigma_*phi*complex_i+d)/(b-rho_*sigma_*phi*complex_i-d);
//
//
//    complex<double> C =  r_*phi*tau*complex_i  + a_/(sigma_*sigma_)*
//								    ( 
//					                   (b-rho_*sigma_*phi*complex_i + d)*tau - 2.0*log( (1.0-g*exp(d*tau))/(1.0-g) )
//				                    );
//
//    complex<double> D = (b-rho_*sigma_*phi*complex_i +d)/(sigma_*sigma_) * (1.0-exp(d*tau))/(1.0-g*exp(d*tau));
//	 
//	complex<double> ff  = exp(C + D*v + x*phi*complex_i);
//
//	return ff;
//}


complex<double> Heston_analytical_pricer :: f(double S, double v,
											  double b, double u,
											  double phi)
{

	double tau = T_;
	double x = log(S);

	complex<double> d1 = sqrt( 
					 (rho_*sigma_*phi*complex_i-b)*(rho_*sigma_*phi*complex_i-b) - sigma_*sigma_*(2*u*phi*complex_i - phi*phi)
				    );


	complex<double> g1 = (b-rho_*sigma_*phi*complex_i + d1)/(b-rho_*sigma_*phi*complex_i-d1);
	//complex<double> g2 = complex_unity/g1;

	////! first formula 
	complex<double> d = d1; 
	complex<double> g = g1;
   
	//! second formula 
	if ( OPTIMIZATION_if_use_CF_second_formula == true)
	{
		d = -d1; 
		g = complex_unity/g1;
	}
   

    complex<double> C =  r_*phi*tau*complex_i  + a_/(sigma_*sigma_)*
								    ( 
					                   (b-rho_*sigma_*phi*complex_i + d)*tau - 2.0*log( (1.0-g*exp(d*tau))/(1.0-g) )
				                    );

    complex<double> D = (b-rho_*sigma_*phi*complex_i  + d)/(sigma_*sigma_) * (1.0-exp(d*tau))/(1.0-g*exp(d*tau));
	 
	complex<double> ff  = exp(C + D*v + x*phi*complex_i);

	return ff;
}


double Heston_analytical_pricer :: integrant( double S, double v, 
				 double b, double u,
				 double phi)
{
	 complex<double> result = exp(-complex_i*phi*log(K_))*f(S,v,b,u,phi)/(complex_i*phi);
	 return result.real();
}


double Heston_analytical_pricer :: Rieman_integerate( double S, double v,
						 double b,   double u, 
						 int num_simulation) 
{
	int discret_size = num_simulation;
	double pace = (integral_interval_end_-integral_interval_begin_)/discret_size;
    vector<double> grill(discret_size+1);

	grill[0] = integral_interval_begin_;
	for(unsigned int k=0; k<grill.size()-1; ++k)
		grill[k+1] = grill[k] + pace;

	double sum = 0;

	for(unsigned int k=0; k<grill.size(); ++k)
	{
		double integrant_temp = integrant(S,  v, b,  u,  grill[k]);
		sum += integrant_temp;
	}
    
    return sum*pace;
}


double Heston_analytical_pricer :: cubic_spline_integral(double S, double v,
							                             double b,   double u, 
							                             int num_simulation)
{
    int discret_size = num_simulation;
	double pace = (integral_interval_end_-integral_interval_begin_)/discret_size;
    vector<double> x(discret_size+1);

	x[0] = integral_interval_begin_;
	for(unsigned int k=0; k<x.size()-1; ++k)
		x[k+1] = x[k] + pace;

	vector<double> y(discret_size+1);
	for(unsigned int i=0; i<y.size(); ++i)
	{
	    y[i] = integrant(S,  v,  b,  u,  x[i]);
	}

	//! construct Interpolator 
	Interpolator interplator(x,y);
    return interplator.calculate_integral_approximation();
}

double Heston_analytical_pricer :: gauss_legendre_integral(  double S, double v,
															 double b,   double u, 
															 int num_simulation)
{
	throw("prolem of include");
	/*
    int discret_size = num_simulation;
	double pace = (integral_interval_end_-integral_interval_begin_)/discret_size;

	//! prepare x, w
    vector<double> x(discret_size+1);
	vector<double> w(discret_size+1);
	double x_min = integral_interval_begin_;
	double x_max = integral_interval_end_;
	gausslegendre(x_min, x_max, x, w); 
   
	//! prepare y
	vector<double> y(x.size());
	for(unsigned int i=0; i<x.size(); ++i)
	{
		y[i] = integrant(S,  v,  b,  u,  x[i]);
	}

	//! GaussLegendre integral
	double result = 0.0;
	for(unsigned int i=0; i<x.size(); ++i)
	{
	    result += y[i]*w[i];
	}
	return result;
	*/
}

double Heston_analytical_pricer :: P( double S, double v, 
								      double b, double u,
									  int num_simulation)
{
	////! Riemann integration
    //double integration = Rieman_integerate(S,v,b,u, num_simulation);
	//! Cubic spline integration
	double integration = cubic_spline_integral(S, v, b, u, num_simulation);
	////! GaussLegendre integration
	//double integration = gauss_legendre_integral(S, v, b, u, num_simulation);

	return 0.5 + pi_inverse*integration;
}


//! input: S_0, v_0, T-t
double Heston_analytical_pricer :: heston_call(double S_t, double v_t, /*double tau, double  K,*/ int num_simulation)
{ 
    double u1 = 0.5;
    double u2 = -0.5;

    double b1 = kappa_ + lambda_ - rho_*sigma_;  // lambda == 0 here !
    double b2 = kappa_ + lambda_;

    double p1 =  P(S_t, v_t, b1, u1, num_simulation);
    double p2 =  P(S_t, v_t, b2, u2, num_simulation);

    return S_t*p1 - K_*exp(-r_*T_)*p2;
}


double Heston_analytical_pricer :: heston_put(double S_t, double v_t, /*double tau, double K,*/ int num_simulation)
{
    return  heston_call(S_t, v_t, num_simulation) - (S_t-K_*exp(-r_*T_));
}


// ******************************************************************************
// 
//							     VECTORIZATION 
//
// *******************************************************************************
//! for compare the PDE result with analytical result, need to evaluate Heston-call-put for different S_t, and v_t
//! so write the vectorized formular

void Heston_analytical_pricer :: f_prepare_v(double b, double u, 
								 		     vector<double>& phi_v, 
								 		     vector<complex<double>>& C_v, // function of b, u, phi :)
										     vector<complex<double>>& D_v) // function of b, u, phi :)
{
	double tau = T_;
	//double x = log(S);

	if(phi_v.size()!= C_v.size() || C_v.size() != D_v.size())
	{
	    throw ("Error in function Heston_analytical_pricer :: f_prepare(...), vectors' size don't match");
	}

	for(unsigned int i=0; i<phi_v.size(); ++i)
	{
		double phi = phi_v[i];
		complex<double> d1 = sqrt( 
						 (rho_*sigma_*phi*complex_i-b)*(rho_*sigma_*phi*complex_i-b) - sigma_*sigma_*(2*u*phi*complex_i - phi*phi)
						);
		complex<double> g1 = (b-rho_*sigma_*phi*complex_i + d1)/(b-rho_*sigma_*phi*complex_i-d1);
		//complex<double> g2 = complex_unity/g1;

		complex<double> d = d1;
		complex<double> g = g1;

		if ( OPTIMIZATION_if_use_CF_second_formula == true)
		{
			d = -d1; 
			g = complex_unity/g1;
		}

		complex<double> C =  r_*phi*tau*complex_i  + a_/(sigma_*sigma_)*
										( 
										   (b-rho_*sigma_*phi*complex_i + d)*tau - 2.0*log( (1.0-g*exp(d*tau))/(1.0-g) )
										);

		complex<double> D = (b-rho_*sigma_*phi*complex_i +d)/(sigma_*sigma_) * (1.0-exp(d*tau))/(1.0-g*exp(d*tau));
		C_v[i] = C;
		D_v[i] = D;
	}
}


void Heston_analytical_pricer :: integrant_v( double S, double v, 
											  vector<double>& phi_v, 
					 						  vector<complex<double>>& C_v, // function of b, u, phi :)
											  vector<complex<double>>& D_v,
											  vector<double>& integrant_vector) // function of b, u, phi :)
{
	 if(S<0.00000001)
		 S = 0.0000001;  // to avoid ln(0.0) error 
     double x = log(S); 

	 for(unsigned int i=0; i<integrant_vector.size(); ++i)
	 {
	     complex<double> f    = exp(C_v[i] + D_v[i]*v + x*phi_v[i]*complex_i);
		 integrant_vector[i] = (exp(-complex_i*phi_v[i]*log(K_))*f/(complex_i*phi_v[i])).real();
	 }
}

void Heston_analytical_pricer :: P_v(    const vector<double>& S,   const vector<double>& v, 
										 double b,   double u, 
										 Matrix& result, 
										 int num_simulation)
{
	int discret_size = num_simulation;
	double pace = (integral_interval_end_-integral_interval_begin_)/discret_size;
    vector<double> x(discret_size+1, 0.0);

	x[0] = integral_interval_begin_;
	for(unsigned int k=0; k<x.size()-1; ++k)
		x[k+1] = x[k] + pace;

	vector<complex<double>> C_v(x.size());
	vector<complex<double>> D_v(x.size());
	
	//! precalculation ...
    f_prepare_v(b, u, x, C_v, D_v);
			

	vector<double> y(x.size(), 0.0); // integrant_v

	for(unsigned int i=0; i<S.size(); ++i)
	{
		cout.precision(2);
		cout << "analytical pricing complet: " << i/(double)S.size() << "%"<< endl;
		for(unsigned int j=0; j<v.size(); ++j)
		{
			integrant_v(S[i], v[j], x, C_v, D_v, y);
			Interpolator interplator(x, y);  // lots of copy coller ---- todo, make it more efficient ???
			double integration = interplator.calculate_integral_approximation();
			result[i][j] = 0.5 + pi_inverse * integration;
		}
	}
}

//! accelerate the calculation: suppose all the Direction_Parameters are fixed and only change S_0, and v_0
Matrix Heston_analytical_pricer :: heston_call_v(vector<double> S_v, vector<double> v_v, int num_simulation)
{
    double u1 = 0.5;
    double u2 = -0.5;

    double b1 = kappa_ + lambda_ - rho_*sigma_;  // lambda == 0 here !
    double b2 = kappa_ + lambda_;

    Matrix p1_m((int)S_v.size(), (int)v_v.size(), 0.0);
	Matrix p2_m((int)S_v.size(), (int)v_v.size(), 0.0);

	P_v(S_v, v_v, b1, u1, p1_m, num_simulation);
	P_v(S_v, v_v, b2, u2, p2_m, num_simulation);
     
	Matrix result_m((int)S_v.size(), (int)v_v.size(), 0.0);

	double coeff = K_*exp(-r_*T_);
	for(unsigned int i=0; i<S_v.size(); ++i)
	{
		for(unsigned int j=0; j<v_v.size(); ++j)
		{
			result_m[i][j] = S_v[i]*p1_m[i][j] - coeff*p2_m[i][j];
		}
	}

	//! don't trust the analytical price for S == 0: 
	for(unsigned int i=0; i<S_v.size(); ++i)
	{
		if(fabs(S_v[i])<0.00000000001)
		{
			for(unsigned int j=0; j<v_v.size(); ++j)
			{
				//double S0 = S_v[i] ;
				//double K = Option_params_K;
				//double r = Heston_params_r;
				//double T = Option_params_T;
				//! use formula: S - K*e^{-rT}
				result_m[i][j] = 0;
			}
		}
	}

	return result_m;
}


//! accelerate the calculation: suppose all the Direction_Parameters are fixed and only change S_0, and v_0
Matrix Heston_analytical_pricer :: heston_put_v(vector<double> S_v, vector<double> v_v, int num_simulation)
{

	Matrix price_m = heston_call_v(S_v, v_v, num_simulation);
    
	for(unsigned int i=0; i<S_v.size(); ++i)
	{
		double call_put_parity_adjust = -(S_v[i]-K_*exp(-r_*T_));
		for(unsigned int j=0; j<v_v.size(); ++j)
		{
			price_m[i][j]  = price_m[i][j] + call_put_parity_adjust;
		}
	}

	return price_m;
}