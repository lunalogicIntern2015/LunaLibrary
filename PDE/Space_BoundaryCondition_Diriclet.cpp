#include "Space_BoundaryCondition_Diriclet.h"

//! initialize_cadre &  get_bc_extremValue & finalize_cadre are the same for Direchlet BC.
double Space_BoundaryCondition_Diriclet :: initialize_cadre(  double t_n, double t_nPlus,
														int index_x, int index_y,
														Matrix& U_0) const
{
	return get_bc_extremValue( t_n, t_nPlus,
								index_x,  index_y,
								U_0);
}

double Space_BoundaryCondition_Diriclet :: finalize_cadre( double t_n, double t_nPlus,
													   int index_x, int index_y,
													   Matrix& U_0) const
{
	return get_bc_extremValue( t_n, t_nPlus,
								index_x,  index_y,
								U_0);
}

//! Adjust the A,B,F,G
//! F and G are zero_column_vector
void Space_BoundaryCondition_Diriclet :: bc_adjust(double t_n, double t_nPlus,
											 int index_x, int index_y,
											 TridiagonalMatrix& A, Matrix& G,
											 TridiagonalMatrix& B, Matrix& F,
											 Matrix& U_0) const
{
	//! Implicit part: G
	double x = discret.get_discret_x(index_x);
	double y = discret.get_discret_y(index_y);

	if(if_leftBoundary == true)          // s_min
	{
		G[0][0] += A.get_left_up_element()*f_diriclet_bc(t_nPlus,x,y);
	}
	else
	{
		G[G.rows()-1][0] += A.get_right_bottom_element()*f_diriclet_bc(t_nPlus,x,y);   // s_max
	}

	//! Explicit part: F 
	double bc_t_n = -99999; 
	if(BC_direction == Direction_Parameters::BC_direction_x)
	{
	    if(if_leftBoundary == true)
		{
		    bc_t_n = U_0[0][index_y];  // not sure if should use U_1 or U_2
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
		throw("Error in function BoundaryCondition_Diriclet :: bc_adjust(...), BC_direction invalide value");
	}

	if(if_leftBoundary == true)
	{
		F[0][0] += B.get_left_up_element()*bc_t_n;           // s_min , bc_L_t_n
	}
	else
	{
		F[F.rows()-1][0] += B.get_right_bottom_element()*bc_t_n;  // s_max,  bc_R_t_n
	}
}


double Space_BoundaryCondition_Diriclet :: get_bc_extremValue(double t_n, double t_nPlus,
														int index_x, int index_y,
														Matrix& U_0) const
{
	double x = discret.get_discret_x(index_x);
	double y = discret.get_discret_y(index_y);

	return f_diriclet_bc(t_nPlus,x,y); 
}

