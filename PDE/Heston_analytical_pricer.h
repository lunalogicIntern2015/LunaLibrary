#pragma once 

//void test_complex();
//double heston_call(double S_t, double v_t, double tau, double  K, int num_simulation);
//double heston_put(double S_t, double v_t, double tau, double K, int num_simulation);
//
#include <complex>  //#include <boost/tr1/complex.hpp>
#include <PDE/Matrix/Matrix.h>
#include <PDE/Params.h>



class Heston_analytical_pricer
{
public:
	 //! Hestion's Params
	 double kappa_;
	 double theta_;
	 double sigma_;
	 double rho_;
	 double lambda_;
	 double r_;
	 double a_;

	  //! Option's Params
	 double K_;
	 double T_;

	 //! Numerical integral;s Params
	 double integral_interval_begin_;
	 double integral_interval_end_;


	 Heston_analytical_pricer();
	 Heston_analytical_pricer(double K);

	 virtual ~ Heston_analytical_pricer(){};

	 complex<double> f( double S, double v,
					    double b, double u,
					    double phi);

	 ////! second formular for characteristic function
	 //complex<double> f2(double S, double v,
		//				  double b, double u,
		//				  double phi);

	 double integrant(double S, double v, 
				 double b, double u,
				 double phi/*, double K*/);

	 double Rieman_integerate(double S, double v,
						 double b,   double u, 
						  int num_simulation) ;

	 double cubic_spline_integral( double S, double v,
							       double b,   double u, 
							       int num_simulation);

	 double gauss_legendre_integral(  double S, double v,
									  double b,   double u, 
									  int num_simulation);

	 double P( double S, double v,
			   double b, double u,
			   int num_simulation);

	 double heston_call(double S_t, double v_t, int num_simulation = Heston_params_num_analytical_fomula_discretization);
	 double heston_put(double S_t, double v_t,  int num_simulation = Heston_params_num_analytical_fomula_discretization);



	 //! VECTORIZATION
	 void  f_prepare_v( double b, double u, 
	 				    vector<double>& phi_v, 
	 				    vector<complex<double>>& C_v, // function of b, u, phi :)
					    vector<complex<double>>& D_v); // function of b, u, phi :)

	 void integrant_v(   double S, double v, 
						  vector<double>& phi_v, 
 						  vector<complex<double>>& C_v, // function of b, u, phi :)
						  vector<complex<double>>& D_v,
						  vector<double>& integrant_vector);

	 void  P_v(  const vector<double>& S,   const vector<double>& v, 
				  double b,   double u, 
				  Matrix& result, 
				  int num_simulation);

	 Matrix heston_call_v(vector<double> S_v, vector<double> v_v, int num_simulation = Heston_params_num_analytical_fomula_discretization);
	 Matrix heston_put_v(vector<double> S_v, vector<double> v_v, int num_simulation = Heston_params_num_analytical_fomula_discretization);
};
