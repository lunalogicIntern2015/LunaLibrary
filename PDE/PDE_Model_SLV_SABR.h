#pragma once

#include "PDE_2D_Model.h"
#include <boost/shared_ptr.hpp>

//! dS = r*S*dt + v*S^{beta}*L*dW1  //when beta = 1, the case of XiuXiu :) 
//! dv = alpha*v*dW2
class PDE_Model_SLV_SABR : public PDE_2D_Model
{
public:
	SABRModel	   sabr_model_;
	Discretization discret_;
	mutable Matrix         L_matrix_; //! to save all the result!  first dimension t, second dimension x!  // this should be updated for each step ! 
	mutable vector<double> L_;        //! for current step, need to chang parameter "double t" --> "int t", if want to put 2 together ...
	std::string    model_discret_type_; 

	//! cache to speed up
	vector<double> x_pow_beta;
	vector<double> x_pow_2_beta;
	Matrix         A_x_helper;

	PDE_Model_SLV_SABR(const SABRModel& sabr_model, const Discretization& discret)
		: sabr_model_(sabr_model),
		discret_(discret),
		L_matrix_(discret.get_sizeDiscret_t(), discret.get_sizeDiscret_x(), 0.0),
		L_(discret.get_sizeDiscret_x(),0.0),
		model_discret_type_("SLV"),
		x_pow_beta(discret.get_sizeDiscret_x(),0.0),
		x_pow_2_beta(discret.get_sizeDiscret_x(),0.0),
		A_x_helper(discret.get_sizeDiscret_x(),discret.get_sizeDiscret_y(),0.0)

	{
		double beta = sabr_model_.get_beta();
		for(unsigned int i=0; i<x_pow_beta.size(); ++i)
		{
			double x = discret_.get_discret_x(i);
			//double y = discret_.get_discret_y( y_index);
			x_pow_beta[i]   = pow(x,beta);
			x_pow_2_beta[i] = x_pow_beta[i]*x_pow_beta[i];
		}

		for( int i=0; i<A_x_helper.rows(); ++i)
		{
			for( int j=0; j<A_x_helper.cols(); ++j)
			{
				double x = discret_.get_discret_x(i);
				double y = discret_.get_discret_y(j);
				A_x_helper[i][j] = 0.5*y*y*x_pow_2_beta[i];
			}
		}

	};

	virtual ~PDE_Model_SLV_SABR(){};

	//! local volativlity function :) 
	double local_vol(double t, double x) const { double r = sabr_model_.get_r();
	return 0.2*x+0.1; //1.0/(x+2.0)*exp(-r*t)*x+0.1; 
	}

	//double A_x(double t, int x_index, int y_index)   const { double x = discret_.get_discret_x( x_index);
	//														 double y = discret_.get_discret_y( y_index);
	//														 double beta = sabr_model_.get_beta();
	//														 return 0.5*y*y*pow(x,2*beta)*L_[x_index]*L_[x_index]; 

	//														 //return 0.5*y*y*x_pow_2_beta[x_index]*L_[x_index]*L_[x_index]; 
	//													   }

	double A_x(double t, int x_index, int y_index)   const { 
		return A_x_helper[x_index][y_index]*L_[x_index]*L_[x_index]; 
	}

	double B_x(double t, int x_index, int y_index)   const { double x = discret_.get_discret_x( x_index);
	//double y = discret_.get_discret_y( y_index);
	double r = sabr_model_.get_r();
	return r*x;
	}

	double C_x(double t, int x_index, int y_index)   const {throw("Error in SABR_Model_discret::C_x(...), not implemented ! ");} 

	double A_y(double t, int x_index, int y_index)   const { //double x = discret_.get_discret_x( x_index);
		double y = discret_.get_discret_y( y_index);
		double alpha = sabr_model_.get_alpha();
		return 0.5*alpha*alpha*y*y;
	}

	double B_y(double t, int x_index, int y_index)   const {return 0;}

	double C_y(double t, int x_index, int y_index)   const {throw("Error in SABR_Model_discret::C_y(...), not implemented ! ");} 

	double F_x_y(double t, int x_index, int y_index) const {  double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double alpha = sabr_model_.get_alpha();
	double beta  = sabr_model_.get_beta();         
	double rho   = sabr_model_.get_rho();
	//if(L_[x_index]!=1.0)
	//{
	// cout << "L = "<< L_[x_index] << endl;
	// throw ("Error L_ error ! ");
	//}
	//return alpha*y*y*pow(x,beta)*rho*L_[x_index];
	return alpha*y*y*x_pow_beta[x_index]*rho*L_[x_index];
	}

	std::string get_model_discret_type() const {return model_discret_type_;}

	//! Matrix D2 size: (matrix_size.first+2, matrix_size.second+2)
	void calculate_L(int t_index,  const Matrix& U_ip) const //! called when calculating from t_index --> t_plus_index = t_index+1
	{
		int t_plus_index = t_index+1;

		double t1 = discret_.get_discret_t(t_index);
		double t2 = discret_.get_discret_t(t_plus_index);
		//! copy slow ... 
		Matrix U(U_ip);
		double integral_nominator   = 0.0; 
		double integral_denominator = 0.0;

		//! integral on v
		for(int i=0; i<discret_.get_sizeDiscret_x(); ++i)
		{
			int x_index = i;
			double x = discret_.get_discret_x(x_index);
			for(int j=0; j<discret_.get_sizeDiscret_y(); ++j)
			{
				double v        = discret_.get_discret_y(j); 
				double density  = U[i][j];
				double interval = 0.0;
				if(j==0) // left-most & right-most
				{
					interval = discret_.get_delta_x(0)/2;
				}
				else if(discret_.get_sizeDiscret_y()-1)
				{
					interval = discret_.get_delta_x(discret_.get_sizeDiscret_y_tilde()-1)/2;
				}
				else
				{
					interval = (discret_.get_delta_x(j-1)+discret_.get_delta_x(j))/2;
				}
				integral_nominator	   += density*interval;
				integral_denominator   += v*v*density*interval;
			}
			if(integral_denominator ==0)
			{
				L_[x_index]                      = 0.0;
				L_matrix_[t_plus_index][x_index] = 0.0;
			}
			else
			{
				L_[x_index]                      = sqrt(pow(local_vol(t2,x),2)*integral_nominator/integral_denominator);
				L_matrix_[t_plus_index][x_index] = L_[x_index]; // backup !
			}
		}
	}
};




//////! dS = v*S^beta*dW1
//////! dv = alpha*v*dW2
////class SLV_Heston_Model_discret : public Model_discret
////{
////public:
////	HestonModel heston_model_;
////
////	SLV_Heston_Model_discret(HestonModel& heston_model):heston_model_(heston_model){};
////
////	virtual ~SLV_Heston_Model_discret(){};
////
////	double L(double t, double x)  {return 0;}
////
////	double A_x(double t, double x, double y)   const {  double l = L(t,x);
////														return 0.5*y*x*x*l*l;}
////
////	double B_x(double t, double x, double y)   const {double r = heston_model_.get_r(); 
////	                                                  return r*x;}
////
////	double C_x(double t, double x, double y)   const {throw("Error in SLV_Heston_Model_discret::C_x(...), not implemented ! ");} 
////
////	double A_y(double t, double x, double y)   const {double sigma = heston_model_.get_sigma();
////													  return 0.5*sigma*sigma*y;}
////
////	double B_y(double t, double x, double y)   const {double kappa = heston_model_.get_kappa();
////													  double theta = heston_model_.get_theta();
////													  return kappa*(theta-y);}
////
////	double C_y(double t, double x, double y)   const {throw("Error in SLV_Heston_Model_discret::C_x(...), not implemented ! ");} 
////
////	double F_x_y(double t, double x, double y) const {double rho = heston_model_.get_rho();
////													  double sigma = heston_model_.get_sigma(); 
////													  return rho*sigma*y*x;}
////};
////
