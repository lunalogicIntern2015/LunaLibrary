#pragma once

#include "PDE_2D_Model.h"
#include <boost/shared_ptr.hpp>
//! dS = v*S^{beta}*dW1
//! dv = alpha*v*dW2
class PDE_Model_SABR : public PDE_2D_Model
{
public:
	SABRModel sabr_model_;
	Discretization discret_;
	std::string    model_discret_type_; 

	PDE_Model_SABR(const SABRModel& sabr_model, const Discretization& discret)
		: sabr_model_(sabr_model),
		discret_(discret),
		model_discret_type_("Normal")
	{};

	virtual ~PDE_Model_SABR(){};

	double A_x(double t, int x_index, int y_index)   const {    double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double beta = sabr_model_.get_beta();
	return 0.5*y*y*pow(x,2*beta);}

	double B_x(double t, int x_index, int y_index)  const { return 0;}

	double C_x(double t, int x_index, int y_index)  const { throw("Error in SABR_Model_discret::C_x(...), not implemented ! ");} 

	double A_y(double t, int x_index, int y_index)  const { double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double alpha = sabr_model_.get_alpha();
	return 0.5*alpha*alpha*y*y;}

	double B_y(double t, int x_index, int y_index)  const { return 0;}

	double C_y(double t, int x_index, int y_index)    const { throw("Error in SABR_Model_discret::C_y(...), not implemented ! ");} 

	double F_x_y(double t, int x_index, int y_index)  const { double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double alpha = sabr_model_.get_alpha();
	double beta  = sabr_model_.get_beta();         
	double rho   = sabr_model_.get_rho();
	return alpha*y*y*pow(x,beta)*rho;}

	std::string get_model_discret_type() const {return model_discret_type_;}
};
