#pragma once

#include "TimeInitial_BoundaryCondition.h"
#include "Space_BoundaryCondition.h"
#include "Space_BoundaryCondition_Diriclet.h"
#include "Space_BoundaryCondition_Neumann.h"
#include "Space_BoundaryCondition_PDE.h"

#include "BoundaryCondition_D2.h"
#include "useful_function.h"
//#include "HestonModel.h"
#include <PDE/HestonModel.h>

#include <iostream>

//! it should be a tempelate class of model!

class BoundaryCondition_Heston_Put : public BoundaryCondition_D2
{
public:
	HestonModel    heston_model_;
	int ffd;  // ??? what's this ???

	//! Inner class
	class BC_Heston_Put_X_L : public Space_BoundaryCondition_Diriclet
	{
		Term_Structure* cts_; 
		PayoffYuan*     payoff_;
	public:
		 BC_Heston_Put_X_L( bool if_leftBoundary, 
							const Direction_Parameters::Direction_Parameters BC_direction,
							const Discretization& discret_ip,
							Term_Structure* cts_ip, 
							PayoffYuan* payoff_ip)						   
			: Space_BoundaryCondition_Diriclet(if_leftBoundary, BC_direction, discret_ip),
			  cts_(cts_ip), 
			  payoff_(payoff_ip)
		 {};

		 virtual ~BC_Heston_Put_X_L(){};

		 double f_diriclet_bc (double t, double x, double y) const
		 {
			  return payoff_->get_strike()*cts_->get_DiscountFactor(t);  //! K*exp(-r*t)
		 }
	};

	class BC_Heston_Put_X_R : public Space_BoundaryCondition_Neumann
	{
	public:
		 BC_Heston_Put_X_R( bool if_leftBoundary,
							const Direction_Parameters::Direction_Parameters BC_direction,
							const Discretization& discret_ip)
			  : Space_BoundaryCondition_Neumann(if_leftBoundary, BC_direction, discret_ip){};
		 virtual ~BC_Heston_Put_X_R(){};

		 double f_neumann_bc (double t, double x, double y) const
		 {
			  return 0;
		 }

	};

	class BC_Heston_Put_Y_L : public Space_BoundaryCondition_PDE
	{
	public:
		//! To have access to the outer class member
		//! C++ not like java don't have access to the outer class attribute automatically
		BoundaryCondition_Heston_Put* outer;

		BC_Heston_Put_Y_L( BoundaryCondition_Heston_Put* outer_ip, //! pointer to its outer class :)
						   bool if_leftBoundary,
						   const Direction_Parameters::Direction_Parameters BC_direction,
						   const Discretization& discret_ip)
			: Space_BoundaryCondition_PDE(if_leftBoundary, BC_direction, discret_ip),
			  outer(outer_ip) {};

		virtual ~BC_Heston_Put_Y_L(){};
		
		double PDE_coeff_U   	 (double t_n, double t_nPlus, double x, double y) const {return -outer->heston_model_.get_r();}													
		double PDE_coeff_t		 (double t_n, double t_nPlus, double x, double y) const {return 1.0;}												
		double PDE_coeff_x_1	 (double t_n, double t_nPlus, double x, double y) const {return x*outer->heston_model_.get_r();}							
		double PDE_coeff_x_2	 (double t_n, double t_nPlus, double x, double y) const {return 0.0;}													
		double PDE_coeff_y_1	 (double t_n, double t_nPlus, double x, double y) const {return outer->heston_model_.get_theta()*outer->heston_model_.get_kappa();}  
		double PDE_coeff_y_2     (double t_n, double t_nPlus, double x, double y) const {return 0.0;}      
		double PDE_coeff_crossing(double t_n, double t_nPlus, double x, double y) const {return 0.0;}
	};

	class BC_Heston_Put_Y_R : public Space_BoundaryCondition_Neumann
	{
	public:
		 BC_Heston_Put_Y_R( bool if_leftBoundary,
							const Direction_Parameters::Direction_Parameters BC_direction,
							const Discretization& discret_ip)
			  : Space_BoundaryCondition_Neumann(if_leftBoundary, BC_direction, discret_ip){};
		 virtual ~BC_Heston_Put_Y_R(){};

		 double f_neumann_bc(double t, double x, double y) const
		 {
			  return 0;
		 }
	};

	//! this function should not be used at all.
	class BC_Heston_Put_I : public TimeInitial_BoundaryCondition
	{
		PayoffYuan* payoff_;
	public:
		 BC_Heston_Put_I( PayoffYuan* payoff_ip)
			  : TimeInitial_BoundaryCondition(){};
		 virtual ~BC_Heston_Put_I(){};

		 double get_inital_bc(double x, double y) const
		 {
			 return max(payoff_->get_strike()-x,0);
		 }
	};

	 //! Constructor
	 BoundaryCondition_Heston_Put(const HestonModel& heston_model,
								  const Discretization& discret,
								  Term_Structure* cts_ip, 
								  PayoffYuan* payoff_ip)
		 : heston_model_(heston_model)
	 {
		 Direction_Parameters::Direction_Parameters direction_x = Direction_Parameters::BC_direction_x;
		 Direction_Parameters::Direction_Parameters direction_y = Direction_Parameters::BC_direction_y;

         bc_X_L_ = Space_BoundaryCondition_PTR(new BC_Heston_Put_X_L(true , direction_x, discret, cts_ip, payoff_ip)); //! Dirichlet
		 bc_X_R_ = Space_BoundaryCondition_PTR(new BC_Heston_Put_X_R(false, direction_x, discret));                    //! Neumann
		 bc_Y_L_ = Space_BoundaryCondition_PTR(new BC_Heston_Put_Y_L(this, true , direction_y, discret));			  //! Diriclet
		 bc_Y_R_ = Space_BoundaryCondition_PTR(new BC_Heston_Put_Y_R(false, direction_y, discret));                    //! Neumann
		 bc_I_   = TimeInitial_BoundaryCondition_PTR(new BC_Heston_Put_I (payoff_ip));
	 }

	 virtual ~BoundaryCondition_Heston_Put(){}
};

