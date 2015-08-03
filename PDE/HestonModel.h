#pragma once

#include "Model.h"

class HestonModel : public Model
{
public:
    double r_;
	double sigma_;
	double kappa_;
	double theta_;
	double rho_;
	double lambda_;

	HestonModel(double r, double sigma, double kappa, double theta, double rho, double lambda = 0)
	{	   r_      = r;
		   sigma_  = sigma;
		   kappa_  = kappa;
		   theta_  = theta;
		   rho_    = rho;
		   lambda_ = lambda;
	}
	~HestonModel(){};
	
	//! To be considered latter 
	double A_s(double t, double s, double v) const {return 0.5*v*s*s;}
	double B_s(double t, double s, double v) const {return r_*s;}
	double C_s(double t, double s, double v) const {return -r_;}

	double A_v(double t, double s, double v) const {return kappa_*(theta_-v) - lambda_;}
	double B_v(double t, double s, double v) const {return 0;}
	double C_v(double t, double s, double v) const {return rho_*sigma_*v*s;}

	//double F_x_y(double t, double s, double v) const {return rho_*sigma_*v*s;}

	//! TODO ... ...
	//! replace all these functions by fucntion "get_Direction_Parameters()" which return a vector!
	double get_r()      const {return r_;}
	double get_sigma()  const {return sigma_;}
	double get_kappa()  const {return kappa_;}
	double get_theta()  const {return theta_;}
	double get_rho()    const {return rho_;}
	double get_lambda() const {return lambda_;}
};