#pragma once 

#include "Space_BoundaryCondition.h"
#include <boost/shared_ptr.hpp>

class Space_BoundaryCondition_Diriclet : public Space_BoundaryCondition
{ 
	 //bool if_left_or_right_;  //! left: ==true, right ==false
public:
     
	 Space_BoundaryCondition_Diriclet( bool            if_leftBoundary_ip,
								 const Direction_Parameters::Direction_Parameters& BC_direction_ip,
								 const Discretization& discret_ip)
		 : Space_BoundaryCondition(if_leftBoundary_ip,
									 BC_direction_ip,
									 discret_ip)
	 {}
     
	 virtual ~Space_BoundaryCondition_Diriclet(){};

	 virtual double f_diriclet_bc(double t, double x, double y) const= 0;

	 double initialize_cadre( double t_n, double t_nPlus,
									int index_x, int index_y,
									Matrix& U_0) const;
	 //! Adjust the A,B,F,G
	 void bc_adjust( double t_n,     double t_nPlus,
					 int    index_x, int index_y, 
					 //double bc_t_n,  // known bc of time t_n
					 TridiagonalMatrix& A, Matrix& G,
					 TridiagonalMatrix& B, Matrix& F,
					 Matrix& U_0) const;

     double get_bc_extremValue(double t_n, double t_nPlus, 
								int index_x,      int index_y,
								Matrix& U_0) const;

	 double finalize_cadre( double t_n, double t_nPlus,
							int index_x, int index_y,
							Matrix& U_0) const;
};


typedef boost::shared_ptr<Space_BoundaryCondition_Diriclet> Space_BoundaryCondition_Diriclet_PTR;
typedef boost::shared_ptr<const Space_BoundaryCondition_Diriclet> Space_BoundaryCondition_Diriclet_CONSTPTR;
