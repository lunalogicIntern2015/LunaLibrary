//#include "Scheme_ADI.h"
#include "Linear_Equation_Solver.h"
#include "Scheme_Yanenko.h"


Scheme_Yanenko :: Scheme_Yanenko(
					const Discretization& discret,
					//Model_discret& model_discret,
					PDE_2D_Model& model_discret,
					BoundaryCondition_D2_CONSTPTR bc
				  )
				  :
				  Scheme_ADI(discret,
							 //model_discret,
							 model_discret,
							 bc),
				  Y         (matrix_size.first+2,matrix_size.second+2,0),
				  U_result  (matrix_size.first+2,matrix_size.second+2,0)
{}

void Scheme_Yanenko :: construct_B_x_tri  (double t1, double t2, int index_y) 
{
	B_x_tri_ = TridiagonalMatrix(matrix_size.first, "id");  
}


void Scheme_Yanenko :: construct_A_x_tri  (double t1, double t2, int index_y) 
{
	construct_L_x_nPlus(t2, index_y);
    double scalar = -(t2-t1); //delta_t_;
    A_x_tri_ = L_x_nPlus_*scalar;  // Not efficient at all :)
	A_x_tri_.addIdMatrix();
}

void Scheme_Yanenko :: initializeScheme_x  (double t1, double t2, int index_y, Matrix& U_0)
{
	construct_B_x_tri(t1, t2, index_y); // use "id", 
	construct_A_x_tri(t1, t2, index_y);
	//! Attension: caluclate_explicit_cross_derivative() should be called before it! 
	initlialize_F_G_x(t1, t2, index_y, U_0);

    //! L_0
	get_cross_derivative_x(H_x_t_n_, index_y);
	H_x_t_n_ *= 0.5*(t2-t1); //delta_t_;
	F_x_t_n_ += H_x_t_n_;
}


void Scheme_Yanenko :: construct_B_y_tri  (double t1, double t2, int index_x) 
{
	B_y_tri_ = TridiagonalMatrix(matrix_size.second, "id");  
}

void Scheme_Yanenko :: construct_A_y_tri  (double t1, double t2, int index_x) 
{
	construct_L_y_nPlus(t2, index_x);

    double scalar = -(t2-t1); //delta_t_;
    A_y_tri_ = L_y_nPlus_*scalar;
	A_y_tri_.addIdMatrix();
}

void Scheme_Yanenko :: initializeScheme_y  (double t1, double t2, int index_x, Matrix& U_0)
{
	construct_B_y_tri(t1, t2, index_x); // use "id"
	construct_A_y_tri(t1, t2, index_x);
	//! Attension: caluclate_explicit_cross_derivative() should be called before it! 
	initlialize_F_G_y(t1, t2, index_x, U_0);

	//! L_0
	get_cross_derivative_y(H_y_t_n_, index_x);
	H_y_t_n_ *= 0.5*(t2-t1); //delta_t_;
	F_y_t_n_ += H_y_t_n_;
}


void Scheme_Yanenko::calculate_one_step_x (int index_t_n, int index_y, Matrix& U_0)
{
	double t_n          = discret_.get_discret_t(index_t_n);
	double t_nPlus      = discret_.get_discret_t(index_t_n+1);

    double t1			= t_n; 
	double t2			= t_nPlus;

	if(index_y==0)  //! on the cadre:  v==v_min  --> set directly the value! 
	{
		//cout << "The first boundary condition must can be set as ExtremValue ! " << endl;
	    for(int index_itr_x=0; index_itr_x<discret_.get_sizeDiscret_x(); ++index_itr_x)
		{
			Y[index_itr_x][index_y] = bc_->bc_Y_L_->initialize_cadre(t1, t2,
						  										   index_itr_x, index_y,//x, y, 
																   U_0);
		}
	}
	else if(index_y==U_0.cols()-1)  //! on the cadre:  v== v_max --> set directly the value! 
	{
		for(int index_itr_x=0; index_itr_x<discret_.get_sizeDiscret_x(); ++index_itr_x)
		{
			Y[index_itr_x][index_y] = bc_->bc_Y_R_->finalize_cadre(t1,t2,
																 index_itr_x, index_y,//x, y,
																 U_0);
		}
	}
	else  // solve 1D PDE 
	{
		//! U^n_x: copy index_y'eme column of U --> U_x
		copy_matrixColumnToVector_withoutFirstAndLastElment(U_0,U_1_x,index_y); // without 1st & last element

		//! Construct A,B,H
		initializeScheme_x(t1, t2, index_y, U_0);

		//! solve linear Equation:  A*U^{n+1/2} + G^{n+1/2} = B*U^{n}+ F^{n} --> U^{n+1/2} other than the boundary value
		Linear_Equation_Solver::resolve(get_A_x_tri_ref(), // A_x_ will be changed! 
									    get_B_x_tri_ref(),
									    F_x_t_n_,
									    G_x_t_nPlus_,
									    U_1_x,		    // input
									    U_2_x,
										U_x_temp);		    // output

		copy_vectorToMatrixColumn_withoutFirstAndLastElment(Y, U_2_x, index_y); //! copy Uplus_x --> index_y'eme column of Uplus

		//! Get the boudnary value of U^{n+1/2}
		int x_left_border_index  = 0;
		int x_right_border_index = discret_.get_sizeDiscret_x()-1;

		//! if it is the Dirichlet BC or Neumann BC, should use U_2, 
		//! PDE BC not sure ...  
		//! Maybe it should be treated before solving the Linear_equation ??? TODO ... 
		Y[0][index_y]            = bc_->bc_X_L_->get_bc_extremValue( t1, t2, x_left_border_index,  index_y, U_0);//, U_2);
		Y[U_0.rows()-1][index_y] = bc_->bc_X_R_->get_bc_extremValue( t1, t2, x_right_border_index, index_y, U_0);//, U_2);
	}
}

void Scheme_Yanenko::calculate_one_step_y (int index_t_n, int index_x, Matrix& U_0)
{

	double t_n          = discret_.get_discret_t(index_t_n);
	double t_nPlus      = discret_.get_discret_t(index_t_n+1);

	double t1			= t_n; 
	double t2			= t_nPlus;


	if(index_x==0)   // x==x_min boundary condition
	{
		for(int index_itr_y=0; index_itr_y<discret_.get_sizeDiscret_y(); ++index_itr_y)
		{
			U_result[index_x][index_itr_y] = bc_->bc_X_L_->initialize_cadre(t1,t2,
																	 index_x,index_itr_y, 
																	 U_0);
		}
	}
	else if(index_x==U_0.rows()-1)  // x=x_max boundary condition
	{
		for(int index_itr_y=0; index_itr_y<discret_.get_sizeDiscret_y(); ++index_itr_y)
		{
			U_result[index_x][index_itr_y] = bc_->bc_X_R_->finalize_cadre(t1,t2,
																	   index_x,index_itr_y, //x, y,
																	   U_0);
		}
	}
	else  // PDE
	{
		//! Construct A,B,F,G,H
		initializeScheme_y(t1, t2, index_x, U_0);

		//! UplusHalf_y: copy index_x'eme row of UplusHalf --> UplusHalf_y  
		copy_matrixRowToVector_withoutFirstAndLastElment(Y, U_1_y, index_x);

		//! solve linear Equation:  A*U^{n+1} + G^{n+1} = B*U^{n+1/2}+ F^{n+1/2} --> U^{n+1} other than the boundary value
		Linear_Equation_Solver::resolve(  get_A_y_tri_ref(), // A_x_ will be changed! 
										  get_B_y_tri_ref(),
										  F_y_t_n_,
										  G_y_t_nPlus_,
										  U_1_y,           // input
										  U_2_y,
										  U_y_temp);		   // output

		copy_vectorToMatrixRow_withoutFirstAndLastElment(U_result, U_2_y, index_x); //! copy U_x --> index_x'eme row of U

		int y_left_border_index  = 0;
		int y_right_border_index = discret_.get_sizeDiscret_y()-1;

		U_result[index_x][0]			  = bc_->bc_Y_L_->get_bc_extremValue(t1, t2, index_x, y_left_border_index,  U_0);
		U_result[index_x][U_0.cols()-1]   = bc_->bc_Y_R_->get_bc_extremValue(t1, t2, index_x, y_right_border_index, U_0);
	}
}

void Scheme_Yanenko :: calculate_one_step(int index_t, Matrix& U_0)
{
	    Y = U_0;
		double t_n	       = discret_.get_discret_t(index_t);         // in fact "t" = tau = T - "t"
		double t_nPlus     = discret_.get_discret_t(index_t+1);

		caluclate_explicit_cross_derivative(t_n, U_0);

		//! total_index_y = 1 to before last
		for(int index_y = 0; index_y<matrix_size.second+2; ++index_y)
		{
			calculate_one_step_x(index_t, index_y, U_0);
		}
		//cout << "U after x directon" << endl;
		//Uplus.print();
		//system("Pause");


		////! ---- OUTPUT 3 ---- 
		//if(DEBUG_printing_result == true)
		//{
		//	std::stringstream ss_i_U_x;
		//	ss_i_U_x << index_t;
		//	ss_i_U_x << "_";
		//	ss_i_U_x << (matrix_size.second+2 -1);
		//	std::string outputFile_U_x = DEBUG_output_path + "X\\UplusHalfe_x_" + ss_i_U_x.str() + ".csv";
		//	print_result(outputFile_U_x, scheme_.discret_.get_discret_x(),scheme_.discret_.get_discret_y(),  Uplus);
		//}


		//! ***********************************************************************
		//! *						   Cross-Derivative for Y
		//! ***********************************************************************
		caluclate_explicit_cross_derivative(t_nPlus, Y);  // ???????  or t_n or t_Plus ???????? 

		////! ---- OUTPUT 4 ---- 
		//if(DEBUG_printing_result == true)
		//{
		//	std::stringstream ss_i_cross_y;
		//	ss_i_cross_y << i;
		//	std::string file_name_crossderivative_y = DEBUG_output_path + "CROSS\\Cross-Derivative_y_" + ss_i_cross_y.str() + ".csv";
		//	print_result(file_name_crossderivative_y, discret.get_discret_x_tilde(),discret.get_discret_y_tilde(),  scheme_.CrossDerivative_);
		//}


		//! ***********************************************************************
		//! *						Y Direction 
		//! ***********************************************************************
		//! Iterate x, and do y-direction calculation: t_nPlusHalf (Uplus)--> t_nPlus (U)

        //std::cout << "Evaluation: T = "<< t_n << "  --> T = " << t_nPlus << "  Direction: y " << std::endl;
		//! index_x from 1 to before last
		//Uplus.print();
		for(int index_x=0; index_x<matrix_size.first+2; ++index_x)
		{
			//! ATTENTION Uplus <--> U exchang role
			calculate_one_step_y(index_t, index_x, U_0);
		}

		U_0 = U_result; // update the result for the last step!
}