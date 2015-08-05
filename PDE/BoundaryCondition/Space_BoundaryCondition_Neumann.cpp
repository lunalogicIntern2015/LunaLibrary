
//-----------------------------------------------------
//! 
//! Bad Bad boundary condition class !!!!! ??????
//!
//-----------------------------------------------------


//! Adjust the A,B,F,G
//! F and G are zero_column_vector

	//! row[0] = upper(i-j=-1), the last element is added and it equals to 0
	//! row[1] = diagonal(i-j=0),
	//! row[2] = lower(i-j=1), the first elemen is added and it equals to 0

#include <PDE/BoundaryCondition/Space_BoundaryCondition_Neumann.h>
#include <boost/shared_ptr.hpp>

//! initialize_cadre  are not implemented, because I don't know how to do that... 
//! get_bc_extremValue & finalize_cadre are the same for Direchlet BC.
double Space_BoundaryCondition_Neumann :: initialize_cadre(double t_n, double t_nPlus,
													 int index_x, int index_y,
													 Matrix& U_0) const
{
    throw ("Error in function BoundaryCondition_Neumann :: initialize_cadre(...), method not implemented, because I don't know how to do it");
}

double Space_BoundaryCondition_Neumann :: finalize_cadre( double t_n, double t_nPlus,
												   int index_x, int index_y,
												   Matrix& U_0) const
{
    return get_bc_extremValue( t_n, t_nPlus,
								index_x,  index_y,
								U_0);
}

void Space_BoundaryCondition_Neumann :: bc_adjust( double t_n, double t_nPlus,
											 int index_x, int index_y,
											 TridiagonalMatrix& A, Matrix& G,
											 TridiagonalMatrix& B, Matrix& F,
											 Matrix& U_0)const 
{
	//! Impliciit part: G 
	double x = discret.get_discret_x(index_x);  // suppose discretization (x,y direction) don't chang with time.
	double y = discret.get_discret_y(index_y);

	if(if_leftBoundary == true)	// s_min
	{
		A.tridiagonalMatrix_[1][0] += A.get_left_up_element();                        // TridiagonalMatrix's left_up
		G[0][0]     = -A.get_left_up_element()*f_neumann_bc(t_nPlus,x,y)*epsilon_BC_neumann;   
	}
	else                        // s_max
	{
		A.tridiagonalMatrix_[1][A.get_sizeMatrix()-1] += A.get_right_bottom_element();               // TridiagonalMatrix's right_bottom
		G[G.rows()-1][0]      = A.get_right_bottom_element()*f_neumann_bc(t_nPlus,x,y)*epsilon_BC_neumann;
	}

	//! Explicit part: F  
	double bc_t_n = -99999; 
	if(BC_direction == Direction_Parameters::BC_direction_x)
	{
	    if(if_leftBoundary == true)
		{
		    bc_t_n = U_0[0][index_y];  // ???? not sure if should use U_1 or U_2
		}
		else
		{
		    bc_t_n = U_0[U_0.rows()-1][index_y];
		}
	}
	else if(BC_direction == Direction_Parameters::BC_direction_y)
	{
		if(if_leftBoundary == true)
		{
		    bc_t_n = U_0[index_x][0];
		}
		else
		{
		    bc_t_n = U_0[index_x][U_0.cols()-1];
		}
	}
	else
	{
		throw("Error in function BoundaryCondition_Neumann :: bc_adjust(...), BC_direction invalide value");
	}
	if(if_leftBoundary == true)
	{
		F[0][0] = B.get_left_up_element()*bc_t_n;           // s_min, bc_L_t_n
	}
	else
	{
		F[F.rows()-1][0] = B.get_right_bottom_element()*bc_t_n;  // s_max, bc_R_t_n
	}

}


double Space_BoundaryCondition_Neumann :: get_bc_extremValue(double t_n,  double t_nPlus,  // called only when the corresponding 1D PDE is solved ...
													   int index_x, int index_y,
													   Matrix& U_0) const
{
	//! Newmann need the implicit boundary condition ...
	double U_2nd     = -99999;            // after 1st element 
	double U_Nmois2  = -99999; // before last element

	if(BC_direction == Direction_Parameters::BC_direction_x)
	{
		//! for x col's 1 and before last element
	     U_2nd     = U_0[1][index_y];            // after 2nd element    // Why use U_2 before ??? 
	     U_Nmois2  = U_0[U_0.rows()-2][index_y]; // before last element  // Why use U_2 bedore ???
	}
	else if(BC_direction == Direction_Parameters::BC_direction_y)
	{
		//! for y row's 1 and before last element
	     U_2nd     = U_0[index_x][1];            // after 2nd element    // Why use U_2 before ??? 
	     U_Nmois2  = U_0[index_x][U_0.cols()-2]; // before last element	 // Why use U_2 before ???
	}
	else
	{
	    throw ("Error in function BoundaryCondition_Neumann :: get_bc_extremValue(...), BC_direction invalide value");
	}

	double x = discret.get_discret_x(index_x);
	double y = discret.get_discret_y(index_y);
	if(if_leftBoundary == true)	// s_min
	{
	     return U_2nd - f_neumann_bc(t_nPlus,x,y)*epsilon_BC_neumann;       // U1!!! U_0 = U_1 - \alpha*\delta
	}
	else                        // s_max
	{
	     return U_Nmois2 + f_neumann_bc(t_nPlus,x,y)*epsilon_BC_neumann;  // U_Nmois2!!! U_1 = U_0 + \alpha*\delta
	}
}

