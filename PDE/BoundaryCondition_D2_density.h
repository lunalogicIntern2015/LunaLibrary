#pragma once

#include "Space_BoundaryCondition_Diriclet.h"
#include "Space_BoundaryCondition_Neumann.h"
#include "Space_BoundaryCondition_PDE.h"
#include "TimeInitial_BoundaryCondition.h"
#include "Space_BoundaryCondition.h"
#include "BoundaryCondition_D2.h"
#include "useful_function.h"
//#include "HestonModel.h"
#include "Model.h"
#include <boost/shared_ptr.hpp>
#include <iostream>

//! it should be a tempelate class of model!

class BoundaryCondition_D2_density : public BoundaryCondition_D2
{
public:
	Model    model_;
	double   density_coeff;

	double get(){return 0;}

	//! Inner class
	class BC_density_X_L : public Space_BoundaryCondition_Diriclet
	{
	public:
		 BC_density_X_L(  bool if_leftBoundary, 
							const Direction_Parameters::Direction_Parameters BC_direction,
							const Discretization& discret_ip)
						
			: Space_BoundaryCondition_Diriclet(if_leftBoundary, BC_direction, discret_ip){};
		 virtual ~BC_density_X_L(){};

		 double f_diriclet_bc (double t, double x, double y) const
		 {
			  //cout << "ldsfjdj ~~~ " << payoff->get_strike() << "  , " << cts->get_DiscountFactor(t) << endl;
			  //getchar();
			  //return payoff->get_strike()*cts->get_DiscountFactor(t);  //! K*exp(-r*t)

			  //! Fwd PDE
			  return 0;
		 }
	};

	class BC_density_X_R : public Space_BoundaryCondition_Diriclet
	{
	public:
		 BC_density_X_R( bool if_leftBoundary,
							const Direction_Parameters::Direction_Parameters BC_direction,
							const Discretization& discret_ip)
			  : Space_BoundaryCondition_Diriclet(if_leftBoundary, BC_direction, discret_ip){};
		 virtual ~BC_density_X_R(){};

		 double f_diriclet_bc (double t, double x, double y) const
		 {
			  return 0;
		 }

	};

	class BC_density_Y_L : public Space_BoundaryCondition_Diriclet
	{
	public:
		//! To have access to the outer class member
		//! C++ not like java don't have access to the outer class attribute automatically
		BoundaryCondition_Heston_Put* outer;

		BC_density_Y_L(  bool if_leftBoundary,
						   const Direction_Parameters::Direction_Parameters BC_direction,
						   const Discretization& discret_ip)
						   : Space_BoundaryCondition_Diriclet( if_leftBoundary, BC_direction, discret_ip){};

		virtual ~BC_density_Y_L(){};
		double f_diriclet_bc (double t, double x, double y) const
		{
			  return 0;
		}
	};

	class BC_density_Y_R : public Space_BoundaryCondition_Diriclet
	{
	public:
		 BC_density_Y_R(bool if_leftBoundary,
							const Direction_Parameters::Direction_Parameters BC_direction,
							const Discretization& discret_ip)
			  : Space_BoundaryCondition_Diriclet( if_leftBoundary, BC_direction, discret_ip){};
		 virtual ~BC_density_Y_R(){};
		 double f_diriclet_bc (double t, double x, double y) const
		 {
			  return 0;
		 }
	};

	//! this function should not be used at all.
	class BC_density_I : public TimeInitial_BoundaryCondition
	{
		Term_Structure* cts_; 
		PayoffYuan* payoff_;

	public:
		double S0_;
		double v0_;
		double init_density_coeff_;
	public:
		 BC_density_I(double S0,
					  double v0,
					  const Discretization& discret,
					  Term_Structure* cts_ip, 
					  PayoffYuan* payoff_ip)
			  : TimeInitial_BoundaryCondition(),
			    cts_(cts_ip),
			    payoff_(payoff_ip),
			    S0_(S0),
				v0_(v0),
				init_density_coeff_(1)
		 {
			 bool flag_x = false;
			 bool flag_y = false;
			 for(int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)  // suppose it is not the first or last grid !
			 {
				 int x_index = i+1;
				 double x = discret.get_discret_x(x_index);
			     if(fabs(x-S0_)<epsilon_compare)
				 {
					 double interval_left  = discret.get_delta_x(i);
					 double interval_right = discret.get_delta_x(i+1);
				     init_density_coeff_  /= (interval_right + interval_left)/2.0;
					 flag_x = true;
					 break;
				 }
			 }

			 for(int i=0; i<discret.get_sizeDiscret_y_tilde(); ++i)
			 {
				 int y_index = i+1;
			 	 double y = discret.get_discret_y(y_index);
			     if(fabs(y-v0_)<epsilon_compare)
				 {
					 double interval_left  = discret.get_delta_y(i);
					 double interval_right = discret.get_delta_y(i+1);
				     init_density_coeff_  /= (interval_right + interval_left)/2.0;
					 flag_y = true;
					 break;
				 }
			 }
			 if(flag_x == false && flag_y == false)
			 {
			     throw ("Error in Boundary condition for Fwd PDE, BC_density_I(), initial value not found in grid!");
			 }
		 };

		 virtual ~BC_density_I(){};

		 double get_inital_bc(double x, double y) const
		 {
			 // return max(payoff->get_strike()-x,0);
             //if(x==S0_ && y == v0_)
			 if(fabs(x-S0_)<epsilon_compare && fabs(y-v0_)<epsilon_compare)
			 {
				 return init_density_coeff_;
			 }
			 else
			 {
			     return 0;
			 }
		 }
	};

	 //! Constructor
	 BoundaryCondition_D2_density(const Model& model,
								  const Discretization& discret,
								  double S0,
								  double v0,
								  Term_Structure* cts_ip, 
								  PayoffYuan* payoff_ip)
		 : model_(model)
	 {
		 Direction_Parameters::Direction_Parameters direction_x = Direction_Parameters::BC_direction_x;
		 Direction_Parameters::Direction_Parameters direction_y = Direction_Parameters::BC_direction_y;

		 bc_X_L_ = Space_BoundaryCondition_PTR(new BC_density_X_L(true,  direction_x, discret)); 
         bc_X_R_ = Space_BoundaryCondition_PTR(new BC_density_X_R(false, direction_x, discret)); //! Neumann
		 bc_Y_L_ = Space_BoundaryCondition_PTR(new BC_density_Y_L(true , direction_y, discret)); //! Diriclet
		 bc_Y_R_ = Space_BoundaryCondition_PTR(new BC_density_Y_R(false, direction_y, discret)); //! Neumann

		 bc_I_   = TimeInitial_BoundaryCondition_PTR(new BC_density_I  (S0, v0, discret, cts_ip, payoff_ip));
	 }

	 virtual ~BoundaryCondition_D2_density(){}
};

typedef boost::shared_ptr<BoundaryCondition_D2_density> BoundaryCondition_D2_density_PTR;
typedef boost::shared_ptr<const BoundaryCondition_D2_density> BoundaryCondition_D2_density_CONSTPTR;
