#pragma once
#include <iostream>
#include <PDE/Term_Structure.h>
#include <PDE/Payoff.h>
#include <PDE/Matrix/TridiagonalMatrix.h>
#include <PDE/Discretization/Discretization.h>
#include <PDE/Matrix/Matrix.h>
#include <PDE/Params.h>
#include <boost/shared_ptr.hpp>


class Space_BoundaryCondition
{
public:
	 // Recuper pointer, neither creat nor really copy them.
 	 //Term_Structure* cts;    // r_t
	 //PayoffYuan* payoff;  

	 bool if_leftBoundary;    // true-->upBoundary(s_min), flase-->downBoundary(s_max)
	 Direction_Parameters::Direction_Parameters BC_direction;
	 Discretization  discret;



     //! construcotr of the abstract class
	 Space_BoundaryCondition(bool if_leftBoundary_ip,
							   const Direction_Parameters::Direction_Parameters& BC_direction_ip,
							   const Discretization& discret_ip)
		 :   if_leftBoundary(if_leftBoundary_ip),
			 BC_direction(BC_direction_ip),
		     discret(discret_ip){};

	 virtual ~Space_BoundaryCondition(){};

	 //! Adjust the A,B,F,G
	 //! F and G are zero_column_vector
	 //! s can be x or y depending on the situation!
	 virtual double initialize_cadre( double t_n, double t_nPlus,
									int index_x, int index_y,
									Matrix& U_0) const = 0;

	 virtual void bc_adjust( double t_n, double t_nPlus,
							 int index_x,  int  index_y, 
							 TridiagonalMatrix& A, Matrix& G, 
							 TridiagonalMatrix& B, Matrix& F,
							 Matrix& U_0) const = 0; 

	 
     virtual double get_bc_extremValue(  double t_n,  double t_nPlus, 
										 int index_x, int index_y,
										 Matrix& U_0) const = 0;

	 virtual double finalize_cadre( double t_n, double t_nPlus,
									int index_x, int index_y,
									Matrix& U_0) const = 0; // for 4 extrem border: y = y_min,y_max;  x = x_min, x_max


	 //! working flow: 
	 //! 1. initialize_cadre. 
	 //! 2. for (int i=1; i<size-1; ++i)
	 //!    {
	 //!        2.1 bc_adjust --> 2.2 solving PDE --> 2.3 get_bc_extremValue
	 //!    }
	 //! 3. finalize_cadre.

	 //! Attention: bc_adjust & get_bc_extremValue should be coherent !

};

typedef boost::shared_ptr<Space_BoundaryCondition> Space_BoundaryCondition_PTR;
typedef boost::shared_ptr<const Space_BoundaryCondition> Space_BoundaryCondition_CONSTPTR;



