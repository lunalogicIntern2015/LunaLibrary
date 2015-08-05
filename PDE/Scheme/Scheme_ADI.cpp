#include <PDE/Scheme/Scheme_ADI.h>
#include <PDE/Solver/Linear_Equation_Solver.h>
#include <PDE/Params.h> //for Fwd PDE

//const bool DEBUG_if_use_half_time = false;  // This is a problem ... !!!!! ??????

//! --------------------------------------------------------------
//! 
//!     ATTENTION: We suppose time discretization is uniform ! 
//! 
//! --------------------------------------------------------------


//!  A*U_2 + G = B*U_1 + F + H, where G,F is the boundary condition, and H is other explicit term.  
Scheme_ADI::Scheme_ADI(
				const Discretization& discret,
				PDE_2D_Model& model_discret,
				BoundaryCondition_D2_CONSTPTR bc
			 )
			 :	     
          matrix_size(pair<int,int>(discret.get_sizeDiscret_x_tilde(),discret.get_sizeDiscret_y_tilde())),
          
		  L_x_n_       (matrix_size.first, "id"),
		  B_x_tri_     (matrix_size.first, "id"),
		  L_x_nPlus_   (matrix_size.first, "id"),
		  A_x_tri_     (matrix_size.first, "id"),

		  F_x_t_n_     (matrix_size.first, 1, 0.0),         
		  G_x_t_nPlus_ (matrix_size.first, 1, 0.0), 
		  H_x_t_n_     (matrix_size.first, 1, 0.0),
 
		  L_y_n_       (matrix_size.second, "id"),
		  B_y_tri_     (matrix_size.second, "id"),
		  L_y_nPlus_   (matrix_size.second, "id"),
		  A_y_tri_     (matrix_size.second, "id"),

		  F_y_t_n_     (matrix_size.second, 1, 0.0),         
		  G_y_t_nPlus_ (matrix_size.second, 1, 0.0), 
		  H_y_t_n_     (matrix_size.second, 1, 0.0),

		  CrossDerivative_(matrix_size.first, matrix_size.second, 0.0),

		  //schemasType_(type),
		  discret_(discret),
		  //model_(model),
		 // model_discret_(model_discret),
		  //model_discret_(model_discret),
		  model_discret_(model_discret),
		  bc_(bc),

		  //! working place
		  U_1_x(matrix_size.first ,1,0),
		  U_1_y(matrix_size.second,1,0),
		  
		  U_2_x(matrix_size.first, 1,0),
		  U_2_y(matrix_size.second,1,0),

		  U_x_temp(matrix_size.first ,1,0),
		  U_y_temp(matrix_size.second,1,0)
{
	  //delta_t_ = discret.get_discret_t_pas();
}


//! ****************************************************************
//!		   STEP_0: Calculate the elemente of Matrix: L_x, L_y
//! ****************************************************************
//! common for 2 direction ...
double Scheme_ADI :: a(double h, double h_moins) const
{
	return -h/(h_moins*(h_moins+h));
}

double Scheme_ADI :: b(double h, double h_moins) const
{
	return (h-h_moins)/(h_moins*h); 
}

double Scheme_ADI :: c(double h, double h_moins) const
{
	return h_moins/(h*(h_moins+h));
}

double Scheme_ADI :: e(double h, double h_moins) const
{
	return 2/(h_moins*(h_moins+h));
}

double Scheme_ADI :: f(double h, double h_moins) const
{
	return -2/(h_moins*h);
}

double Scheme_ADI :: g(double h, double h_moins) const
{
	return 2/(h*(h_moins+h));
}



//! ---- ---- Derivative along s ---- ----
//! ----FD
//! D_1 
double Scheme_ADI :: a_x(double t, int x_index) const
{
	double h_moins = discret_.get_delta_x(x_index-1);
	double h = discret_.get_delta_x(x_index); 
	return a(h,h_moins);
}

double Scheme_ADI :: b_x(double t, int x_index) const
{
	double h_moins = discret_.get_delta_x(x_index-1);
	double h = discret_.get_delta_x(x_index); 
	return b(h,h_moins);
}

double Scheme_ADI :: c_x(double t, int x_index) const
{
	double h_moins = discret_.get_delta_x(x_index-1);
	double h = discret_.get_delta_x(x_index);
	return c(h,h_moins);
}
//! D_2
double Scheme_ADI :: e_x(double t, int x_index) const
{
	double h_moins = discret_.get_delta_x(x_index-1);
	double h = discret_.get_delta_x(x_index); 
	return e(h,h_moins);
}
double Scheme_ADI :: f_x(double t, int x_index) const
{
	double h_moins = discret_.get_delta_x(x_index-1);
	double h = discret_.get_delta_x(x_index); 
	return f(h,h_moins);
}

double Scheme_ADI :: g_x(double t, int x_index) const
{
	double h_moins = discret_.get_delta_x(x_index-1);
	double h = discret_.get_delta_x(x_index); 
	return g(h,h_moins);
}
//! ---- L_x
//double Scheme_ADI :: l0_x(int t_index, int x_index, int y_index) const
double Scheme_ADI :: l0_x(double t, int x_index, int y_index) const
{
	//double t = discret_.get_discret_t(t_index);
	if(if_fwd_PDE == false)
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_x(t,x,y)            *e_x(t,x_index) + model_discret_.B_x(t,x,y)              *a_x(t,x_index);
		return model_discret_.A_x(t,x_index,y_index)*e_x(t,x_index) + model_discret_.B_x(t,x_index,y_index)*a_x(t,x_index);
		
	}
	else
	{
		//double x = discret_.get_discret_x(x_index-1);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_x(t,x,y)            *e_x(t,x_index) - model_discret_.B_x(t,x,y)              *a_x(t,x_index);
        int x_index_fwd = x_index -1;
		int y_index_fwd = y_index;
		return model_discret_.A_x(t,x_index_fwd,y_index_fwd)*e_x(t,x_index) - model_discret_.B_x(t,x_index_fwd,y_index_fwd)*a_x(t,x_index);
	}
}

double Scheme_ADI :: l1_x(double t, int x_index, int y_index) const
{
	if(if_fwd_PDE == false)
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_x(t,x,y)            *f_x(t,x_index) + model_discret_.B_x(t,x,y)              *b_x(t,x_index) + model_discret_.C_x(t,x,y);
		return model_discret_.A_x(t,x_index,y_index)*f_x(t,x_index) + model_discret_.B_x(t,x_index,y_index)*b_x(t,x_index) + model_discret_.C_x(t,x_index,y_index);
	}
	else
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_x(t,x,y)			*f_x(t,x_index) - model_discret_.B_x(t,x,y)			   *b_x(t,x_index);
		int x_index_fwd = x_index;
		int y_index_fwd = y_index;
		return model_discret_.A_x(t,x_index_fwd,y_index_fwd)*f_x(t,x_index) - model_discret_.B_x(t,x_index_fwd,y_index_fwd)*b_x(t,x_index);
	}
}

double Scheme_ADI :: l2_x(double t, int x_index, int y_index) const
{
	if(if_fwd_PDE == false)
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_x(t,x,y)			*g_x(t,x_index) + model_discret_.B_x(t,x,y)			   *c_x(t,x_index);
		return model_discret_.A_x(t,x_index,y_index)*g_x(t,x_index) + model_discret_.B_x(t,x_index,y_index)*c_x(t,x_index);
	}
	else
	{
		//double x = discret_.get_discret_x(x_index+1);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_x(t,x,y)			*g_x(t,x_index) - model_discret_.B_x(t,x,y)			   *c_x(t,x_index);
		int x_index_fwd = x_index+1;
		int y_index_fwd = y_index;
		double xxxx = model_discret_.A_x(t,x_index_fwd,y_index_fwd);
		//double yyyy = g_x(t,x_index);
		//double zzzz = model_discret_.B_x(t,x_index_fwd,y_index_fwd);
		//double cccc = c_x(t,x_index);
		return model_discret_.A_x(t,x_index_fwd,y_index_fwd)*g_x(t,x_index) - model_discret_.B_x(t,x_index_fwd,y_index_fwd)*c_x(t,x_index);
	}
}


//! ---- ---- Derivative along y ---- ----
//! ----FD
//! D_1 
double Scheme_ADI :: a_y(double t, int y_index) const
{
	double h_moins = discret_.get_delta_y(y_index-1);
	double h = discret_.get_delta_y(y_index); 
	return a(h,h_moins);
}

double Scheme_ADI :: b_y(double t, int y_index) const
{
	double h_moins = discret_.get_delta_y(y_index-1);
	double h = discret_.get_delta_y(y_index); 
	return b(h,h_moins);
}

double Scheme_ADI :: c_y(double t, int y_index) const
{
	double h_moins = discret_.get_delta_y(y_index-1);
	double h = discret_.get_delta_y(y_index); 
	return c(h,h_moins);
}
//! D_2
double Scheme_ADI :: e_y(double t, int y_index) const
{
	double h_moins = discret_.get_delta_y(y_index-1);
	double h = discret_.get_delta_y(y_index); 
	return e(h,h_moins);
}
double Scheme_ADI :: f_y(double t, int y_index) const
{
	double h_moins = discret_.get_delta_y(y_index-1);
	double h = discret_.get_delta_y(y_index); 
	return f(h,h_moins);
}

double Scheme_ADI :: g_y(double t, int y_index) const
{
	double h_moins = discret_.get_delta_y(y_index-1);
	double h = discret_.get_delta_y(y_index); 
	return g(h,h_moins);
}

//! ---- L_y
double Scheme_ADI :: l0_y(double t, int x_index, int y_index) const
{
	if(if_fwd_PDE == false)
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_y(t,x,y)			*e_y(t,y_index) + model_discret_.B_y(t,x,y)			   *a_y(t,y_index);
		return model_discret_.A_y(t,x_index,y_index)*e_y(t,y_index) + model_discret_.B_y(t,x_index,y_index)*a_y(t,y_index);
	}
	else
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index-1);
		//return model_discret_.A_y(t,x,y)			*e_y(t,y_index) - model_discret_.B_y(t,x,y)			   *a_y(t,y_index);
		int x_index_fwd = x_index;
		int y_index_fwd = y_index-1;
		return model_discret_.A_y(t,x_index_fwd,y_index_fwd)*e_y(t,y_index) - model_discret_.B_y(t,x_index_fwd,y_index_fwd)*a_y(t,y_index);
	}
}

double Scheme_ADI :: l1_y(double t, int x_index, int y_index) const
{
	if(if_fwd_PDE == false)
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index); 	
		//return model_discret_.A_y(t,x,y)			*f_y(t,y_index) + model_discret_.B_y(t,x,y)			   *b_y(t,y_index) + model_discret_.C_y(t,x,y);
		return model_discret_.A_y(t,x_index,y_index)*f_y(t,y_index) + model_discret_.B_y(t,x_index,y_index)*b_y(t,y_index) + model_discret_.C_y(t,x_index,y_index);
	}
	else
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index); 	
		//return model_discret_.A_y(t,x,y)			*f_y(t,y_index) - model_discret_.B_y(t,x,y)			  *b_y(t,y_index);
		int x_index_fwd = x_index;
		int y_index_fwd = y_index;
		return model_discret_.A_y(t,x_index_fwd,y_index_fwd)*f_y(t,y_index) - model_discret_.B_y(t,x_index_fwd,y_index_fwd)*b_y(t,y_index);
	}
}

double Scheme_ADI :: l2_y(double t, int x_index, int y_index) const
{
	if(if_fwd_PDE == false)
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index);
		//return model_discret_.A_y(t,x,y)			*g_y(t,y_index) + model_discret_.B_y(t,x,y)			   *c_y(t,y_index);
		return model_discret_.A_y(t,x_index,y_index)*g_y(t,y_index) + model_discret_.B_y(t,x_index,y_index)*c_y(t,y_index);
	}
	else
	{
		//double x = discret_.get_discret_x(x_index);
		//double y = discret_.get_discret_y(y_index+1);
		//return model_discret_.A_y(t,x,y)			*g_y(t,y_index) - model_discret_.B_y(t,x,y)			   *c_y(t,y_index);
		int x_index_fwd = x_index;
		int y_index_fwd = y_index+1;
		return model_discret_.A_y(t,x_index_fwd,y_index_fwd)*g_y(t,y_index) - model_discret_.B_y(t,x_index_fwd,y_index_fwd)*c_y(t,y_index);
	}
}


//! ****************************************************************
//!		    STEP_1: Explicite Cross Derivative
//! ****************************************************************
void Scheme_ADI :: caluclate_explicit_cross_derivative(double t, Matrix& U_t_n) 
//! U_n's size              = (matrix_size.first+2, matrix_size.second+2) 
//! CrossDerivative_'s size = (matrix_size.first, matrix_size.second)
{
	if(if_fwd_PDE == false)
	{
		for(int i_tilde=0; i_tilde<matrix_size.first; ++i_tilde)
		{
			 for(int j_tilde=0; j_tilde<matrix_size.second; ++j_tilde)
			 {
				 int i = i_tilde + 1;
				 int j = j_tilde + 1;

				 CrossDerivative_[i_tilde][j_tilde]  =   U_t_n[i+1][j+1] 
													   + U_t_n[i-1][j-1]
													   - U_t_n[i-1][j+1]
													   - U_t_n[i+1][j-1];

				 double h_x_moins   = discret_.get_delta_x(i_tilde);
				 double h_x         = discret_.get_delta_x(i_tilde+1);
				 double h_y_moins   = discret_.get_delta_y(j_tilde);
				 double h_y         = discret_.get_delta_y(j_tilde+1);  
				 double denominator = h_x*h_y + h_x_moins*h_y_moins + h_x*h_y_moins + h_x_moins*h_y;

				 CrossDerivative_[i_tilde][j_tilde] /=  denominator;

				 //double x = discret_.get_discret_x(i);
				 //double y = discret_.get_discret_y(j);
				 //CrossDerivative_[i_tilde][j_tilde] *=  model_discret_.F_x_y(t,x,y);

				 int x_index = i;
				 int y_index = j;
				 CrossDerivative_[i_tilde][j_tilde] *=  model_discret_.F_x_y(t,x_index,y_index);
			 }
		}   
		//cout << "CrossDerivative_" << endl;
		//CrossDerivative_.print();
	}
	else
	{
		for(int i_tilde=0; i_tilde<matrix_size.first; ++i_tilde)
		{
			 for(int j_tilde=0; j_tilde<matrix_size.second; ++j_tilde)
			 {
				 int i = i_tilde + 1;
				 int j = j_tilde + 1;
     //            
				 //double x_plus            = discret_.get_discret_x(i+1);
				 //double x_moins           = discret_.get_discret_x(i-1);
				 //double y_plus            = discret_.get_discret_y(j+1);
				 //double y_moins           = discret_.get_discret_y(j-1);

				 int x_index_plus  = i+1;
				 int x_index_moins = i-1;
				 int y_index_plus  = j+1;
				 int y_index_moins = j-1;

				 //double coeff_plus_plus   = model_discret_.F_x_y(t,x_plus, y_plus );
				 //double coeff_plus_moins  = model_discret_.F_x_y(t,x_plus, y_moins);
				 //double coeff_moins_plus  = model_discret_.F_x_y(t,x_moins,y_plus );
				 //double coeff_moins_moins = model_discret_.F_x_y(t,x_moins,y_moins);


				 double coeff_plus_plus   = model_discret_.F_x_y(t,x_index_plus, y_index_plus );
				 double coeff_plus_moins  = model_discret_.F_x_y(t,x_index_plus, y_index_moins);
				 double coeff_moins_plus  = model_discret_.F_x_y(t,x_index_moins,y_index_plus );
				 double coeff_moins_moins = model_discret_.F_x_y(t,x_index_moins,y_index_moins);


				 CrossDerivative_[i_tilde][j_tilde]  =   coeff_plus_plus*  U_t_n[i+1][j+1] 
													   + coeff_moins_moins*U_t_n[i-1][j-1]
													   - coeff_moins_plus* U_t_n[i-1][j+1]
													   - coeff_plus_moins* U_t_n[i+1][j-1];

				 double h_x_moins   = discret_.get_delta_x(i_tilde);
				 double h_x         = discret_.get_delta_x(i_tilde+1);
				 double h_y_moins   = discret_.get_delta_y(j_tilde);
				 double h_y         = discret_.get_delta_y(j_tilde+1);  
				 double denominator = h_x*h_y + h_x_moins*h_y_moins + h_x*h_y_moins + h_x_moins*h_y;

				 CrossDerivative_[i_tilde][j_tilde] /=  denominator;

				 ////double t = discret_.get_discret_t(t_index);
				 //double x = discret_.get_discret_x(i);
				 //double y = discret_.get_discret_y(j);
				 //CrossDerivative_[i_tilde][j_tilde] *=  model_discret_.F_x_y(t,x,y);
			 }
		}   
	}
}
//! Crossing derivative getting a col:
void Scheme_ADI :: get_cross_derivative_x(Matrix& cross_derivative_x, unsigned int index_y) 
{
    //! We don't stock the boudnary element for the cross-derivative 
    int tilde_index_y = index_y - 1;

    //!  TODO: should not test matrix size here but int he matrix class! 
    if(cross_derivative_x.rows() != matrix_size.first || cross_derivative_x.cols() != 1)
    {
	    cout << "Error in function get_cross_derivative_col(...), matrix size dismatch! " << endl;
    }
    else if (tilde_index_y < 0 || tilde_index_y > matrix_size.second-1)
    {
        cout << "Error in function get_cross_derivative_col(...), total_index_y takes no valid value ! "<< endl;
    }
    else
    {
        for( int i=0; i<matrix_size.first; ++i)
	    {
	        cross_derivative_x[i][0] = CrossDerivative_[i][tilde_index_y];
	    } 
    }
}

//! Crossing derivative getting a row:
void Scheme_ADI :: get_cross_derivative_y(Matrix& cross_derivative_y, unsigned int index_x) 
{
    //! We don't stock the boudnary element for the cross-derivative 
    int tilde_index_x = index_x - 1;

    if(cross_derivative_y.rows() != matrix_size.second || cross_derivative_y.cols() != 1)
    {
	    cout << "Error in function get_cross_derivative_row(...), matrix size dismatch! " << endl;
    }
    else if(tilde_index_x<0 || tilde_index_x  > matrix_size.first-1)
    {
        cout << "Error in function get_cross_derivative_col(...), total_index_x takes no valid value ! "<< endl;
    }
    else
    {
       for( int j=0; j<matrix_size.second; ++j)
	   {
	     cross_derivative_y[j][0] = CrossDerivative_[tilde_index_x][j];
	   } 
    }
}

//! ****************************************************************
//!			 STEP_2: Direction x (y fix)
//! ****************************************************************
//! Direction x 
void Scheme_ADI :: construct_L_x_n    (double t1, int y_index)  // t_index = t_n_index
{
	for(int i_tilde=0; i_tilde<matrix_size.first; i_tilde++)
	{
		int x_index = i_tilde+1;
		L_x_n_.tridiagonalMatrix_[0][i_tilde] = l2_x(t1, x_index, y_index);
		L_x_n_.tridiagonalMatrix_[1][i_tilde] = l1_x(t1, x_index, y_index);
		L_x_n_.tridiagonalMatrix_[2][i_tilde] = l0_x(t1, x_index, y_index);
	}
}

void Scheme_ADI :: construct_L_x_nPlus(double t2, int y_index)  // t_index = t_nPlus_index
{
	//double y = discret_.get_discret_y(y_index);
	for(int i_tilde=0; i_tilde<matrix_size.first; i_tilde++)
	{
		//double x = discret_.get_discret_x_tilde(i_tilde);
		int x_index = i_tilde+1;
		//A_x_tri_.tridiagonalMatrix_[0][i_tilde] = l2_x(t2, x_index, y_index);
		//A_x_tri_.tridiagonalMatrix_[1][i_tilde] = l1_x(t2, x_index, y_index);
		//A_x_tri_.tridiagonalMatrix_[2][i_tilde] = l0_x(t2, x_index, y_index);

		L_x_nPlus_.tridiagonalMatrix_[0][i_tilde] = l2_x(t2, x_index, y_index);
		L_x_nPlus_.tridiagonalMatrix_[1][i_tilde] = l1_x(t2, x_index, y_index);
		L_x_nPlus_.tridiagonalMatrix_[2][i_tilde] = l0_x(t2, x_index, y_index);
	}
}

void Scheme_ADI :: initlialize_F_G_x  (double t1, double t2, int index_y, Matrix& U_0) 
{
    //! ---- TODO: not efficient at all ...
	//! Initialize F,G to zero vector  
	init_to_zero_colMatrix(F_x_t_n_);
	init_to_zero_colMatrix(G_x_t_nPlus_);

	int x_left_border_index  = 0;
	int x_right_border_index = discret_.get_sizeDiscret_x()-1;

	bc_->bc_X_L_->bc_adjust(
		t1, t2,
        x_left_border_index, index_y, //x_left_border,y,
		A_x_tri_, G_x_t_nPlus_,
		B_x_tri_, F_x_t_n_, 
		U_0);//, U_2);

	bc_->bc_X_R_->bc_adjust(
		t1, t2, 
		x_right_border_index, index_y, //x_rihgt_border,	y,		
		A_x_tri_, G_x_t_nPlus_,
		B_x_tri_, F_x_t_n_,
		U_0);//, U_2);

	//get_cross_derivative_x(H_x_t_n_, index_y);
	//H_x_t_n_ *= 0.5*delta_t_;
	//F_x_t_n_ += H_x_t_n_;
}


//void Scheme_ADI :: initializeScheme_x  (double t1, double t2, int index_y, Matrix& U_1, Matrix& U_2)
//{
//	construct_B_x_tri(t1, index_y);
//	construct_A_x_tri(t2, index_y);
//	//! Attension: caluclate_explicit_cross_derivative() should be called before it! 
//	initlialize_F_G_x(t1, t2, index_y, U_1, U_2);
//}

//void Scheme_ADI :: calculate_one_step_x (int index_t_n, bool if_first_half_t, int index_y,
//										 Matrix& U_0, Matrix& U_1, Matrix& U_2, int scheme_x_itr) 
//										 //U_0 is for calculating the boundary condition (set: initial, final, F,G)
//										 //U_1 --> U_2 is for calculating the PDE: (set: A,B and calculate: A*U_2+G = B*U_2+F)
//{
//	double t_n          = discret_.get_discret_t(index_t_n);
//	double t_nPlus      = discret_.get_discret_t(index_t_n+1);
//	//double t_nPlusHalf  = (t_n + t_nPlus)/2.0;
//    double t1			= t_n; 
//	double t2			= t_nPlus;
//	//if(DEBUG_if_use_half_time == true)
//	//{
//	//	if(if_first_half_t == true)
//	//	{
//	//	    t1 = t_n;
//	//		t2 = t_nPlusHalf;
//	//	}
//	//	else
//	//	{
//	//		throw ("Error  ---- haha ---- ");
//	//	    t1 = t_nPlusHalf;
//	//		t2 = t_nPlus;
//	//	}
//	//}
//
//	if(index_y==0)  //! on the cadre:  v==v_min  --> set directly the value! 
//	{
//		//cout << "The first boundary condition must can be set as ExtremValue ! " << endl;
//	    for(int index_itr_x=0; index_itr_x<discret_.get_sizeDiscret_x(); ++index_itr_x)
//		{
//			U_2[index_itr_x][index_y] = bc_.bc_Y_L->initialize_cadre(t1, t2,
//						  											 index_itr_x, index_y,//x, y, 
//																	 U_0);//, U_2);
//		}
//
//  //      //! Output for checking the bug ...
//		//string file_name = "C:\\v_min\\hello.csv";
//		//vector<double> x_vector_value(discret_.get_discret_x());
//		//vector<double> y_vector_value(1,-999999); 
//		//Matrix output_boundary(x_vector_value.size(),y_vector_value.size());
//		//for(int i=0; i<x_vector_value.size(); ++i)
//		//{
//		//	output_boundary[i][0] = U_2[i][0];
//		//}
//		//print_result(file_name , x_vector_value, y_vector_value, output_boundary);
//		////system("Pause");
//
//	}
//	else if(index_y==U_2.cols()-1)  //! on the cadre:  v== v_max --> set directly the value! 
//	{
//		for(int index_itr_x=0; index_itr_x<discret_.get_sizeDiscret_x(); ++index_itr_x)
//		{
//			U_2[index_itr_x][index_y] = bc_.bc_Y_R->finalize_cadre(t1,t2,
//																   index_itr_x, index_y,//x, y,
//																   U_0);//, U_2);
//		}
//	}
//	else  // solve 1D PDE 
//	{
//		//! U^n_x: copy index_y'eme column of U --> U_x
//		copy_matrixColumnToVector_withoutFirstAndLastElment(U_1,U_1_x,index_y); // without 1st & last element
//
//		//! Construct A,B,H
//		initializeScheme_x(t1, t2, index_y, U_0);//, U_1, U_2, scheme_x_itr);  // initialize matrix: A,B,H, need to do : F,G
//
//			//cout << " ---- new scheme matrix A " << endl;
//			//A_x_tri_.print();
//			//cout << " ---- new scheme matrix G " << endl;
//			//G_x_t_nPlus_.print();
//			//cout << " ---- new scheme matrix B " << endl;
//			//B_x_tri_.print();
//			//cout << " ---- new scheme matrix F " << endl;
//			//F_x_t_n_.print();
//			//cout << " ---- new scheme matrix H " << endl;
//			//H_x_t_n_.print();
//			//cout << " ---- new scheme matrix U_1_x " << endl;
//			//U_1_x.print();
//
//
//		//! solve linear Equation:  A*U^{n+1/2} + G^{n+1/2} = B*U^{n}+ F^{n} --> U^{n+1/2} other than the boundary value
//		Linear_Equation_Solver::resolve(get_A_x_tri_ref(), // A_x_ will be changed! 
//									    get_B_x_tri_ref(),
//									    F_x_t_n_,
//									    G_x_t_nPlus_,
//									    U_1_x,		    // input
//									    U_2_x);		    // output
//		//cout << " scheme matrix U_2_x " << endl;
//		//U_2_x.print();
//		//getchar();
//		copy_vectorToMatrixColumn_withoutFirstAndLastElment(U_2, U_2_x, index_y); //! copy Uplus_x --> index_y'eme column of Uplus
//
//		//! Get the boudnary value of U^{n+1/2}
//		int x_left_border_index  = 0;
//		int x_right_border_index = discret_.get_sizeDiscret_x()-1;
//
//		//! if it is the Dirichlet BC or Neumann BC, should use U_2, 
//		//! PDE BC not sure ...  
//		//! Maybe it should be treated before solving the Linear_equation ??? TODO ... 
//		U_2[0][index_y]            = bc_.bc_X_L->get_bc_extremValue( t1, t2, x_left_border_index,  index_y, U_0);//, U_2);
//		U_2[U_2.rows()-1][index_y] = bc_.bc_X_R->get_bc_extremValue( t1, t2, x_right_border_index, index_y, U_0);//, U_2);
//	}
//}
//

//! Direction y
void Scheme_ADI :: construct_L_y_n    (double t1, int x_index)  // t_index = t_n_index
{
	//L_y_n_.tridiagonalMatrix_.print();
	for(int j_tilde=0; j_tilde<matrix_size.second; j_tilde++)
	{
		int y_index = j_tilde+1;
		//if(j_tilde==0 || j_tilde ==1 || j_tilde==9)
		//{
		//	cout << "j_tilde = " << j_tilde << endl;
		//    cout << " l2_y = " << l2_y(t1, x_index, y_index) << endl; 
		//	cout << " l1_y = " << l1_y(t1, x_index, y_index) << endl; 
		//	cout << " l0_y = " << l0_y(t1, x_index, y_index) << endl; 
		//}
		L_y_n_.tridiagonalMatrix_[0][j_tilde] = l2_y(t1, x_index, y_index);
		L_y_n_.tridiagonalMatrix_[1][j_tilde] = l1_y(t1, x_index, y_index);
		L_y_n_.tridiagonalMatrix_[2][j_tilde] = l0_y(t1, x_index, y_index);

	}
	//throw ("Error in function Scheme_ADI :: construct_L_y_n (...), ADI don't use this function, B_x_tri is initialized as id matrix in constructor of scheme");
}

void Scheme_ADI :: construct_L_y_nPlus(double t2, int x_index)   // t_index = t_nPlus_index
{
	//double x = discret_.get_discret_y(index_x);
	for(int j_tilde=0; j_tilde<matrix_size.second; j_tilde++)
	{
		//double y = discret_.get_discret_y_tilde(j_tilde);

		//int y_index = j_tilde+1;
		//A_y_tri_.tridiagonalMatrix_[0][j_tilde] = l2_y(t2, x_index, y_index);
		//A_y_tri_.tridiagonalMatrix_[1][j_tilde] = l1_y(t2, x_index, y_index);
		//A_y_tri_.tridiagonalMatrix_[2][j_tilde] = l0_y(t2, x_index, y_index);
		int y_index = j_tilde+1;
		L_y_nPlus_.tridiagonalMatrix_[0][j_tilde] = l2_y(t2, x_index, y_index);
		L_y_nPlus_.tridiagonalMatrix_[1][j_tilde] = l1_y(t2, x_index, y_index);
		L_y_nPlus_.tridiagonalMatrix_[2][j_tilde] = l0_y(t2, x_index, y_index);
	}
}


void Scheme_ADI :: initlialize_F_G_y  (double t1, double t2,
	                                   int index_x, Matrix& U_0) 
{
    //! ---- TODO: not efficient at all ...
	//! Initialize F,G to zero vector  
	init_to_zero_colMatrix(F_y_t_n_);
	init_to_zero_colMatrix(G_y_t_nPlus_);

	int y_left_border_index  = 0;
	int y_right_border_index = discret_.get_sizeDiscret_y()-1;

	bc_->bc_Y_L_->bc_adjust(
		t1, t2,
		index_x,y_left_border_index,
		A_y_tri_, G_y_t_nPlus_,
		B_y_tri_, F_y_t_n_, 
		U_0);//, U_2);

	bc_->bc_Y_R_->bc_adjust(
		t1, t2, 
		index_x,y_right_border_index,
		A_y_tri_, G_y_t_nPlus_,
		B_y_tri_, F_y_t_n_,
		U_0);//, U_2);

	//get_cross_derivative_y(H_y_t_n_, index_x);
	//H_y_t_n_ *= 0.5*delta_t_;
	//F_y_t_n_ += H_y_t_n_;
}

