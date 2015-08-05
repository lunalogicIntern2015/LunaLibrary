//! Adjust the A,B,F,G
//! F and G are zero_column_vector

	//! row[0] = upper(i-j=-1), the last element is added and it equals to 0
	//! row[1] = diagonal(i-j=0),
	//! row[2] = lower(i-j=1), the first elemen is added and it equals to 0

#include <PDE/BoundaryCondition/Space_BoundaryCondition_PDE.h>
#include <PDE/Heston_analytical_pricer.h>

//! initialize_cadre &  get_bc_extremValue & finalize_cadre are the same for Direchlet BC.
double Space_BoundaryCondition_PDE :: initialize_cadre(double t_n, double t_nPlus,
												 int index_x, int index_y,
												 Matrix& U_1) const
{
	return get_bc_extremValue(  t_n, t_nPlus,
								index_x,  index_y,
								U_1);
}

double Space_BoundaryCondition_PDE :: finalize_cadre(double t_n, double t_nPlus,
												 int index_x, int index_y,
												 Matrix& U_1) const
{
	return get_bc_extremValue(  t_n, t_nPlus,
								index_x,  index_y,
								U_1);
}

//! for the moment I suppose: there is non second order derivative: for simplicity :)
//! when it is called we don't tach the cadre ! 
void Space_BoundaryCondition_PDE :: bc_adjust( double t_n, double t_nPlus,
										 int index_x, int index_y,
										 TridiagonalMatrix& A, Matrix& G,
										 TridiagonalMatrix& B, Matrix& F,
										 Matrix& U_0) const
{
	//! implicit: 
	double extream_value = get_bc_extremValue( t_n,  t_nPlus,
											   index_x,  index_y,
 											   U_0);

  
	if(if_leftBoundary == true)
	{
		G[0][0]  += A.get_left_up_element() * extream_value;	 // but why not x ???
	}
	else 
	{
	    G[G.rows()-1][0]  += A.get_left_up_element() * extream_value;	 // but why not x ???
	}


	//! Explicit: F
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

//void BoundaryCondition_PDE :: bc_adjust( double t_n, double t_nPlus,
//										 double x, double y,
//										 double bc_t_n, // known bc of time t_n
//										 TridiagonalMatrix& A, Matrix& G,
//										 TridiagonalMatrix& B, Matrix& F)  
//{
//
//	if(if_leftBoundary == true)	// s_min
//	{
//
//		A.tridiagonalMatrix_[1][0] += A.get_left_up_element();                        // TridiagonalMatrix's left_up
//		//G[0][0]     = -A.get_left_up_element()*f_neumann_bc(t_nPlus,x,y)*epsilon_BC_neumann;   
//	}
//	else                        // s_max
//	{
//		A.tridiagonalMatrix_[1][A.get_sizeMatrix()-1] += A.get_right_bottom_element();               // TridiagonalMatrix's right_bottom
//		//G[G.rows()-1][0]      = A.get_right_bottom_element()*f_neumann_bc(t_nPlus,x,y)*epsilon_BC_neumann;
//	}
//
//	if(if_leftBoundary == true)
//	{
//		//F[0][0] = B.get_left_up_element()*bc_t_n;           // s_min, bc_L_t_n
//		if(B.get_left_up_element()*bc_t_n  !=0)
//		{
//		    throw ("Error in BoundaryCondition_Diriclet :: bc_adjust(...), it should be zero ! ");
//		}
//	}
//	else
//	{
//		//F[F.rows()-1][0] = B.get_right_bottom_element()*bc_t_n;  // s_max, bc_R_t_n
//		if(B.get_right_bottom_element()*bc_t_n  !=0)
//		{
//		    throw ("Error in BoundaryCondition_Diriclet :: bc_adjust(...), it should be zero ! ");
//		}
//	}
//}



//! when it is called we may tach the cadre: so need to distinguish the case ... 

			////! BC: PDE 
			//if(i_x_whenIndexyEquals0<matrix_size.first+2-1)
			//{
			//	HestonModelYuan* h_model = dynamic_cast<HestonModelYuan*>(schemas_->model_);
			//	double BB = h_model->r_*discret_->get_discret_x(i_x_whenIndexyEquals0);
			//	double CC = -h_model->r_;
			//	double GG = h_model->kappa_*h_model->theta_;
			//	double lambda1 = BB* schemas_->delta_t_/schemas_->h_s_;
			//	double lambda2 = GG* schemas_->delta_t_/schemas_->h_v_;
			//	Uplus[i_x_whenIndexyEquals0][0] = (1 - lambda1- lambda2 + CC*schemas_->delta_t_)*U[i_x_whenIndexyEquals0][0]
			//									  + lambda1 * U[i_x_whenIndexyEquals0+1][0]
			//									  + lambda2 * U[i_x_whenIndexyEquals0][1];
			//}
			//else
			//{
			//	Uplus[i_x_whenIndexyEquals0][0] = 0;
			//}

double Space_BoundaryCondition_PDE :: get_bc_extremValue( double t_n,  double t_nPlus,
												    int index_x, int index_y,
												    Matrix& U_1) const
{
	//////! PDE need the explicit one
	////double U_2nd     = -99999;            // after 1st element 
	////double U_Nmois2  = -99999; // before last element

	////if(BC_direction == Direction_Parameters::BC_direction_x)
	////{
	////	//! for x col's 1 and before last element
	////     U_2nd     = U_1[1][index_y];            // after 1st element 
	////     U_Nmois2  = U_1[U_2.rows()-2][index_y]; // before last element
	////}
	////else if(BC_direction == Direction_Parameters::BC_direction_y)
	////{
	////	//! for y row's 1 and before last element
	////     U_2nd     = U_1[index_x][1];            // after 1st element 
	////     U_Nmois2  = U_1[index_x][U_2.cols()-2]; // before last element	
	////}
	////else
	////{
	////    throw ("Error in function BoundaryCondition_Neumann :: get_bc_extremValue(...), BC_direction invalide value");
	////}

	//There is a bug for this boundary condition ... 
	double x = discret.get_discret_x(index_x);
	double y = discret.get_discret_y(index_y);

	//! (t_n,t_nPlus) = (t_n,t_nPlusHalf) or (t_nPlusHalf,t_nPlus)
    double coeff_U        = PDE_coeff_U       (t_n,t_nPlus,x,y);
	double coeff_t        = PDE_coeff_t       (t_n,t_nPlus,x,y);
	double coeff_x_1      = PDE_coeff_x_1     (t_n,t_nPlus,x,y);
	double coeff_x_2      = PDE_coeff_x_2     (t_n,t_nPlus,x,y);
	double coeff_y_1      = PDE_coeff_y_1     (t_n,t_nPlus,x,y);
	double coeff_y_2      = PDE_coeff_y_2     (t_n,t_nPlus,x,y);
	double coeff_crossing = PDE_coeff_crossing(t_n,t_nPlus,x,y);

	if(coeff_x_2 != 0 || coeff_y_2 != 0)
	{
	    throw("Error in function BoundaryCondition_PDE :: bc_adjust(...): don't implement BC for the case second derivative !=0 ");
	}

	double U_element              = U_1[index_x][index_y]; 

	double U_element_x_bump       = -99999;
	double U_element_y_bump       = -99999;
	double U_element_x_bump_sign  = -99999;
	double U_element_y_bump_sign  = -99999;

	double temp_delta_BC_t = t_nPlus - t_n; 
	double temp_delta_BC_x = -999999;
	double temp_delta_BC_y = -999999;

	//!  default method is: x[i+1] - x[i] & y[i+1] - y[i]
	if(index_x == U_1.rows()-1)
	{
	    U_element_x_bump      = U_1[index_x-1][index_y]; 
		U_element_x_bump_sign = -1;
		temp_delta_BC_x       = discret.get_delta_x(index_x-1);

	}
	else
	{
	    U_element_x_bump      = U_1[index_x+1][index_y]; 
		U_element_x_bump_sign = 1;
		temp_delta_BC_x       = discret.get_delta_x(index_x);
	}

	if(index_y == U_1.cols()-1)
	{
	    U_element_y_bump      = U_1[index_x][index_y-1]; 
		U_element_y_bump_sign = -1;
		temp_delta_BC_y       = discret.get_delta_y(index_y-1);
	}
	else 
	{
		U_element_y_bump      = U_1[index_x][index_y+1]; 
		U_element_y_bump_sign = 1;		   
		temp_delta_BC_y       = discret.get_delta_y(index_y);
	}

	//! only works for uniform discretization on space ! 
	//if(if_leftBoundary == true)  
	//{
	//	temp_delta_BC_x = delta_BC_x_left;
	//	temp_delta_BC_y = delta_BC_y_left;
	//}
	//else if(if_leftBoundary == false)  
	//{ 
	//	temp_delta_BC_x = delta_BC_x_right;
	//	temp_delta_BC_y = delta_BC_y_right;
	//}
	//else
	//{
	//    throw("Error in function BoundaryCondition_PDE :: get_bc_extremValue(...), invalide value if_leftBoundary");
	//}

    double part_x     = coeff_x_1*(U_element_x_bump-U_element)/temp_delta_BC_x*U_element_x_bump_sign;
	double part_y     = coeff_y_1*(U_element_y_bump-U_element)/temp_delta_BC_y*U_element_y_bump_sign;
	double part_U     = coeff_U*U_element;
	double part_x_y_U = part_x + part_y + part_U;

	double u_corner_value = U_element + temp_delta_BC_t/coeff_t*part_x_y_U;

	return u_corner_value;	
}
