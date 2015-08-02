#pragma once
#include "TimeInitial_BoundaryCondition.h"
#include "Space_BoundaryCondition.h"
#include "Space_BoundaryCondition_Neumann.h"
#include <boost/shared_ptr.hpp>

class BoundaryCondition_D2
{
public:
	 //! t = 0
	 TimeInitial_BoundaryCondition_PTR bc_I_;
     
	 Space_BoundaryCondition_PTR bc_X_L_;  // x_min
	 Space_BoundaryCondition_PTR bc_X_R_;  // x_max
	 Space_BoundaryCondition_PTR bc_Y_L_;  // y_min
	 Space_BoundaryCondition_PTR bc_Y_R_;  // y_max

	 //! default constructor
	 BoundaryCondition_D2(){}

	 BoundaryCondition_D2(TimeInitial_BoundaryCondition_PTR bc_I,
						  Space_BoundaryCondition_PTR bc_X_L,
						  Space_BoundaryCondition_PTR bc_X_R,
						  Space_BoundaryCondition_PTR bc_Y_L,
						  Space_BoundaryCondition_PTR bc_Y_R)
						  :bc_I_(bc_I), 
						  bc_X_L_(bc_X_L),
						  bc_X_R_(bc_X_R),
						  bc_Y_L_(bc_Y_L),
						  bc_Y_R_(bc_Y_R){}


	 virtual ~BoundaryCondition_D2(){};
};

typedef boost::shared_ptr<BoundaryCondition_D2> BoundaryCondition_D2_PTR;
typedef boost::shared_ptr<const BoundaryCondition_D2> BoundaryCondition_D2_CONSTPTR;

//----------------------------------------------------------------------------------
//
//								FACTORY
//
//----------------------------------------------------------------------------------
class BoundaryCondition_D2_Factory
{
public:
	static BoundaryCondition_D2_PTR create_BC_D2_space1stDerivative0(TimeInitial_BoundaryCondition_PTR bc_I,
																	 const Discretization& discret_ip)
	{
		Direction_Parameters::Direction_Parameters direction_x = Direction_Parameters::BC_direction_x;
		Direction_Parameters::Direction_Parameters direction_y = Direction_Parameters::BC_direction_y;

		Space_BoundaryCondition_PTR bc_X_L = Space_BoundaryCondition_PTR(new Space_BoundaryCondition_Neumann_1stDerivative0(true, direction_x, discret_ip));
		Space_BoundaryCondition_PTR bc_X_R = Space_BoundaryCondition_PTR(new Space_BoundaryCondition_Neumann_1stDerivative0(false, direction_x, discret_ip));
		Space_BoundaryCondition_PTR bc_Y_L = Space_BoundaryCondition_PTR(new Space_BoundaryCondition_Neumann_1stDerivative0(true, direction_y, discret_ip));
		Space_BoundaryCondition_PTR bc_Y_R = Space_BoundaryCondition_PTR(new Space_BoundaryCondition_Neumann_1stDerivative0(false, direction_y, discret_ip));
		
		return BoundaryCondition_D2_PTR(new BoundaryCondition_D2(bc_I, bc_X_L, bc_X_R, bc_Y_L, bc_Y_R));
	}
};