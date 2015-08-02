#pragma once 

//! TODO: Not a good idea to define: "epsilon_BC_neumann = 0.0001" directly!

#include "Space_BoundaryCondition.h"
#include <boost/shared_ptr.hpp>

class Space_BoundaryCondition_PDE : public Space_BoundaryCondition
{ 
public:
	//! constructor & virtual destructor
	Space_BoundaryCondition_PDE(bool if_leftBoundary_ip,
						  const Direction_Parameters::Direction_Parameters& BC_direction_ip,
						  const Discretization& discret_ip) 
	 :  Space_BoundaryCondition(if_leftBoundary_ip,
								  BC_direction_ip,
								  discret_ip)

	{};

	virtual ~Space_BoundaryCondition_PDE(){};

	//double f_PDE_bc(double t_n, double t_nPlus, double x, double y);
	virtual double PDE_coeff_U       (double t_n, double t_nPlus, double x, double y) const = 0;     //  coefficient of U
	virtual double PDE_coeff_t       (double t_n, double t_nPlus, double x, double y) const = 0;     //  coefficient of 1st derivative of t
	virtual double PDE_coeff_x_1     (double t_n, double t_nPlus, double x, double y) const = 0;     //  coefficient of 1st derivative of s
	virtual double PDE_coeff_x_2     (double t_n, double t_nPlus, double x, double y) const = 0;     //  coefficient of 2nd derivative of s
	virtual double PDE_coeff_y_1     (double t_n, double t_nPlus, double x, double y) const = 0;     //  coefficient of 1st derivative of v
    virtual double PDE_coeff_y_2     (double t_n, double t_nPlus, double x, double y) const = 0;     //  coefficient of 2nd derivative of v
	virtual double PDE_coeff_crossing(double t_n, double t_nPlus, double x, double y) const = 0;     //  coefficient of crossing derivative

	double initialize_cadre( double t_n, double t_nPlus,
									int index_x, int index_y,
									Matrix& U_0) const;
	//! Adjust the A,B,F,G
	void bc_adjust(  double t_n, double t_nPlus,
					 int index_x,  int index_y,
					 //double bc_t_n,  // known bc of time t_n
					 TridiagonalMatrix& A, Matrix& G,
					 TridiagonalMatrix& B, Matrix& F, 
					 Matrix& U_0) const;

	double get_bc_extremValue( double t_n,  double t_nPlus,
							   int index_x, int index_y,
							   Matrix& U_0) const;

	double finalize_cadre( double t_n, double t_nPlus,
									int index_x, int index_y,
									Matrix& U_0) const;
};

typedef boost::shared_ptr<Space_BoundaryCondition_PDE> Space_BoundaryCondition_PDE_PTR;
typedef boost::shared_ptr<const Space_BoundaryCondition_PDE> Space_BoundaryCondition_PDE_CONSTPTR;
