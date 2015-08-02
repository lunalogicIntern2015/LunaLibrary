#pragma once

#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <PDE/PDE_2D_Model.h>
#include <boost/shared_ptr.hpp>

//! Cheyette model hiearchy is not well done, so for the moment only treat the cheyetteDD model

//! suppose underlying markovian states: (x_t,y_t).

class PDE_CheyetteDD_Model: public PDE_2D_Model
{
public:
	CheyetteDD_Model  cheyette_model_;
	Discretization    discret_;

	PDE_CheyetteDD_Model(const CheyetteDD_Model& cheyette_Model, const Discretization& discret)
		: cheyette_model_(cheyette_Model),
		  discret_(discret)
	{};

	virtual ~PDE_CheyetteDD_Model(){};

	double A_x(double t, int x_index, int y_index) const 
	{
		double x = discret_.get_discret_x( x_index);
		double y = discret_.get_discret_y( y_index);

		double diffusion_x = cheyette_model_.diffusion_x(t,x,y);
		return 0.5*diffusion_x*diffusion_x;
	}

	double B_x(double t, int x_index, int y_index)   const
	{
		double x = discret_.get_discret_x( x_index);
		double y = discret_.get_discret_y( y_index);

		double drift_x = cheyette_model_.drift_x_Q(t,x,y);
		return drift_x;
	}

	double C_x(double t, int x_index, int y_index)   const 
	{
		double x = discret_.get_discret_x( x_index);
		double r = cheyette_model_.r_t(t,x); 
		return -r; 
	} 

	double A_y(double t, int x_index, int y_index)   const  // diffusion_y = 0
	{
		return 0.0;
	}

	double B_y(double t, int x_index, int y_index)   const 
	{  
		double x = discret_.get_discret_x( x_index);
		double y = discret_.get_discret_y( y_index);
		double drift_y = cheyette_model_.drift_y(t,x,y);

		return drift_y;
	}

	double C_y(double t, int x_index, int y_index)   const 
	{
		return 0.0;  // YY C_y can be changed to -0.5*r, when C_x = -0.5*r
	} 

	double F_x_y(double t, int x_index, int y_index) const   // because: diffusion_y = 0
	{
		return 0.0;
	}
};

typedef boost::shared_ptr<PDE_CheyetteDD_Model> PDE_CheyetteDD_Model_PTR;
typedef boost::shared_ptr<const PDE_CheyetteDD_Model> PDE_CheyetteDD_Model_CONSTPTR;