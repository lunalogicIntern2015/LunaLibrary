#pragma once 

//! Attention: this boundary condition only works for invariant grid on time :)

#include <PDE/BoundaryCondition/Space_BoundaryCondition.h>
#include <boost/shared_ptr.hpp>
class Space_BoundaryCondition_Neumann : public Space_BoundaryCondition
{ 
public:
	 double epsilon_BC_neumann;
	 vector<double> discret_vector;

	 //! constructor & virtual destructor
	 Space_BoundaryCondition_Neumann(bool if_leftBoundary_ip,
							   const Direction_Parameters::Direction_Parameters& BC_direction_ip,
							   const Discretization& discret_ip)
							   // Note : if use "Discretization& discret", I need to specify the direction x or y, so use directly vector<double>
		 :  Space_BoundaryCondition(if_leftBoundary_ip,
									  BC_direction_ip,
									  discret_ip)
	 {
		 if(BC_direction_ip == Direction_Parameters::BC_direction_x)
		 {
			 discret_vector = discret.get_discret_x();
		 }
		 else if(BC_direction_ip == Direction_Parameters::BC_direction_y)
		 {
             discret_vector = discret.get_discret_y();
		 }
		 else
		 {
		     throw ("Error in constructor of BoundaryCondition_Neumann, invalude Direction_Parameters type ");
		 }

	     if(if_leftBoundary == true)  // s_min
		 {
			 double x1 = discret_vector[0];
			 double x2 = discret_vector[1];
		     epsilon_BC_neumann = x2 - x1;
		 }
		 else // s_max
		 {
			 double x1 = discret_vector[discret_vector.size()-2];
			 double x2 = discret_vector[discret_vector.size()-1];
		     epsilon_BC_neumann = x2 - x1;
		 }
	 };
	 virtual ~Space_BoundaryCondition_Neumann(){};

	 virtual double f_neumann_bc(double t, double x, double y) const = 0;

	 double initialize_cadre( double t_n, double t_nPlus,
									int index_x, int index_y,
									Matrix& U_0) const;
	 //! Adjust the A,B,F,G
	 void bc_adjust( double t_n,  double t_nPlus,
					 int index_x, int index_y,
					 //double bc_t_n,  // known bc of time t_n
					 TridiagonalMatrix& A, Matrix& G,
					 TridiagonalMatrix& B, Matrix& F,
					 Matrix& U_0) const;

	 double get_bc_extremValue(double t_n,  double t_nPlus,
							   int index_x, int index_y,
							   Matrix& U_0) const; 

 	 double finalize_cadre( double t_n, double t_nPlus,
									int index_x, int index_y,
									Matrix& U_0) const;
};
typedef boost::shared_ptr<Space_BoundaryCondition_Neumann> Space_BoundaryCondition_Neumann_PTR;
typedef boost::shared_ptr<const Space_BoundaryCondition_Neumann> Space_BoundaryCondition_Neumann_CONSTPTR;


class Space_BoundaryCondition_Neumann_1stDerivative0 : public Space_BoundaryCondition_Neumann
{
public:
	Space_BoundaryCondition_Neumann_1stDerivative0(bool if_leftBoundary_ip,
		const Direction_Parameters::Direction_Parameters& BC_direction_ip,
		const Discretization& discret_ip):
	Space_BoundaryCondition_Neumann(if_leftBoundary_ip, BC_direction_ip, discret_ip)
	{}

	virtual ~Space_BoundaryCondition_Neumann_1stDerivative0(){}

	double f_neumann_bc(double t, double x, double y) const {return 0;}
};
typedef boost::shared_ptr<Space_BoundaryCondition_Neumann_1stDerivative0> Space_BoundaryCondition_Neumann_1stDerivative0_PTR;
typedef boost::shared_ptr<const Space_BoundaryCondition_Neumann_1stDerivative0> Space_BoundaryCondition_Neumann_1stDerivative0_CONSTPTR;


