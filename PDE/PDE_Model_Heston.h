#pragma once

#include <PDE/HestonModel.h>
#include "PDE_2D_Model.h"
#include <boost/shared_ptr.hpp>

class PDE_Model_Heston: public PDE_2D_Model
{
public:
	HestonModel    heston_model_;
	Discretization discret_;
	std::string    model_discret_type_; 

	PDE_Model_Heston(const HestonModel& heston_model, const Discretization& discret)
		: heston_model_(heston_model),
		discret_(discret),
		model_discret_type_("Normal")
	{};

	virtual ~PDE_Model_Heston(){};

	double A_x(double t, int x_index, int y_index)   const { double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	return 0.5*y*x*x;}

	double B_x(double t, int x_index, int y_index)   const { double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double r = heston_model_.get_r(); 
	return r*x;}

	double C_x(double t, int x_index, int y_index)   const { double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double r = heston_model_.get_r(); 
	return -r/2.0;} 

	double A_y(double t, int x_index, int y_index)   const { double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double sigma = heston_model_.get_sigma();
	return 0.5*sigma*sigma*y;}

	double B_y(double t, int x_index, int y_index)   const {    double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double kappa = heston_model_.get_kappa();
	double theta = heston_model_.get_theta();
	return kappa*(theta-y);}

	double C_y(double t, int x_index, int y_index)   const {    double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double r = heston_model_.get_r(); 
	return -r/2.0;} 

	double F_x_y(double t, int x_index, int y_index) const { double x = discret_.get_discret_x( x_index);
	double y = discret_.get_discret_y( y_index);
	double rho = heston_model_.get_rho();
	double sigma = heston_model_.get_sigma(); 
	return rho*sigma*y*x;}

	std::string get_model_discret_type() const {return model_discret_type_;}
};

typedef boost::shared_ptr<PDE_Model_Heston> PDE_Model_Heston_PTR;
typedef boost::shared_ptr<const PDE_Model_Heston> PDE_Model_Heston_CONSTPTR;