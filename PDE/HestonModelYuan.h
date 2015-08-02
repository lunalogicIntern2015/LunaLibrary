//#pragma once
//#include <math.h>
//#include "Model_D2.h"
//
////! Heston Model
////  dS_t/S_t = rdt + sqrt(v_t)dW^1
////  dv_t     = kappa*(theta-v_t)dt + vol*sqrt(v_t)dW^2
////  d(W_t^1,W_t^2) = rho dt
//
//class HestonModelYuan : public Model_D2
//{
//public:
//    double r_;
//	double sigma_;
//	double kappa_;
//	double theta_;
//	double rho_;
//
//	HestonModelYuan(double r, double sigma, double kappa, double theta, double rho)
//	{	   r_     = r;
//		   sigma_ = sigma;
//		   kappa_ = kappa;
//		   theta_ = theta;
//		   rho_   = rho;
//	}
//	~HestonModelYuan(){};
//	
//	//! To be considered latter 
//	double get_drift_1(double t, double x, double y){return r_*x;}
//	double get_vol_1(double t, double x, double y)  {return sqrt(y)*x;}
//	double get_drift_2(double t, double x, double y){return kappa_*(theta_-y);}
//	double get_vol_2(double t, double x, double y)  {return sigma_*sqrt(y);}
//
//	//! TODO ... ...
//	//! replace all these functions by fucntion "get_parameters()" which return a vector!
//	double get_r()    {return r_;}
//	double get_sigma()  {return sigma_;}
//	double get_kappa(){return kappa_;}
//	double get_theta(){return theta_;}
//	double get_rho()  {return rho_;}
//};
//
