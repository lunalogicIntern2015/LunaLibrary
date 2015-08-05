//#include <PDE/Scheme/Scheme_ADI.h>
#include "Scheme_CS.h"
#include <PDE/Solver/Linear_Equation_Solver.h>
#include <PDE/useful_function.h>
#include <sstream>

using namespace std;
//! Scheme: 0,1,2
//! Scheme_0: what ever direction: Y_0 = (I + L_0 + L_1)*U_0
//! Scheme_1: direction_x: (1-theta*L_1)*Y_1 = Y_0 - theta*L_1*U_0 (at this step: U_1 = Y_0, U_2 = Y_1)
//! Scheme_2: direction_y: (1-theta*L_2)*Y_2 = Y_1 - theta*L_2*U_0 (at this step: U_1 = Y_1, U_2 = Y_2) 

const bool if_print_out = false;

Scheme_CS :: Scheme_CS(
					const Discretization& discret,
					PDE_2D_Model& model_discret,
					BoundaryCondition_D2_CONSTPTR bc,
					double theta, 
					double lambda
				  )
				  :
				  Scheme_ADI( discret,
							  model_discret,
							  bc),
				  scheme_itr_size_(3),
				  theta_(theta),
				  lambda_(lambda),

				  Y_0       (matrix_size.first+2,matrix_size.second+2,0),
				  Y_1       (matrix_size.first+2,matrix_size.second+2,0),
				  Y_2       (matrix_size.first+2,matrix_size.second+2,0)
				  {}

//! ----------------------------------------------------
//! 
//!                   Douglas part
//! 
//! ----------------------------------------------------
//! scheme_0
void Scheme_CS :: initializeScheme_0_x(double t1, double t2, int index_y, Matrix& U_0) 
{
	//B_x_tri_
	construct_L_x_n(t1, index_y);
	double scalar = t2-t1; //delta_t_; //YYYY _change
	B_x_tri_ = L_x_n_*scalar;
	A_x_tri_ = TridiagonalMatrix(matrix_size.first, "id");
	
	initlialize_F_G_x(t1, t2, index_y, U_0);
}

void Scheme_CS :: initializeScheme_0_y(double t1, double t2, int index_x, Matrix& U_0) 
{
	//B_y_tri_
	construct_L_y_n(t1, index_x);
	double scalar =  t2-t1; //delta_t_; //YYYY _change
	B_y_tri_ = L_y_n_*scalar;
	A_y_tri_ = TridiagonalMatrix(matrix_size.second, "id");
	initlialize_F_G_y(t1, t2, index_x, U_0);
}

void Scheme_CS:: calculate_one_step_0 (int index_t_n, Matrix& U_0) 
// Total explicit: Y_0 = (I+ L_0 + L_1 + L_2)U_0
// Don't care about its boundary condition !
{
	double t1  = discret_.get_discret_t(index_t_n);
	double t2  = discret_.get_discret_t(index_t_n+1);

	copy_smallMatrix_to_bigMatrix_withOutBoundary(CrossDerivative_,Y_0);
	Y_0  *=  (t2-t1); //delta_t_; //YYYY _change
	Y_0  += U_0;

	//! Don't care about the boudnary condition ! 
	//! x direction for: L_x_n_
	for(int index_y = 1; index_y<matrix_size.second+2-1; ++index_y)
	{
		//! U^n_x: copy index_y'eme column of U --> U_x
		copy_matrixColumnToVector_withoutFirstAndLastElment(U_0,U_1_x,index_y); // without 1st & last element

		//! Construct A,B,F,G,H
		initializeScheme_0_x(t1, t2, index_y, U_0);

		//! solve linear Equation:  A*U^{n+1/2} + G^{n+1/2} = B*U^{n}+ F^{n} --> U^{n+1/2} other than the boundary value
		Linear_Equation_Solver::resolve(get_A_x_tri_ref(), // A_x_ will be changed! 
										get_B_x_tri_ref(),
										F_x_t_n_,
										G_x_t_nPlus_,
										U_1_x,		    // input
										U_2_x,
										U_x_temp);		    // output

		add_vectorToMatrixColumn_withoutFirstAndLastElment(Y_0, U_2_x, index_y); //! copy Uplus_x --> index_y'eme column of Uplus
	    
		//! Don't care about the boundary condition !
	}

		//! x direction for: L_x_n_
	for(int index_x = 1; index_x<matrix_size.first+2-1; ++index_x)
	{

		//! UplusHalf_y: copy index_x'eme row of UplusHalf --> UplusHalf_y  
		copy_matrixRowToVector_withoutFirstAndLastElment(U_0, U_1_y, index_x);

		//! Construct A,B,F,G,H
		initializeScheme_0_y(t1, t2, index_x, U_0);//, U_1, U_2, scheme_y_itr);	

		//cout << " get_B_y_tri_ref; " << endl;  get_B_y_tri_ref().print();

		//! solve linear Equation:  A*U^{n+1} + G^{n+1} = B*U^{n+1/2}+ F^{n+1/2} --> U^{n+1} other than the boundary value
		Linear_Equation_Solver::resolve(  get_A_y_tri_ref(), // A_x_ will be changed! 
										  get_B_y_tri_ref(),
										  F_y_t_n_,
										  G_y_t_nPlus_,
										  U_1_y,           // input
										  U_2_y,
										  U_y_temp);		   // output

		//cout << " get_A_y_tri_ref; " << endl;  get_A_y_tri_ref().print();

		//get_A_y_tri_ref().print()
		//cout << " U_1_y; " << endl;  U_1_y.print();
		//cout << " U_2_y.print(); " << endl;  U_2_y.print();
		add_vectorToMatrixRow_withoutFirstAndLastElment(Y_0, U_2_y, index_x); //! copy U_x --> index_x'eme row of U
	   
		//! Don't care about the boundary condition !
	}
}


//! scheme_1 direction x :
void Scheme_CS :: initializeScheme_1(double t1, double t2, int index_y, Matrix& U_0)
{
	double scalar = -theta_* (t2-t1); //*delta_t_; //YYYY _changedelta_t_;

	//! B_x_tri_
	construct_L_x_n(t1, index_y);
	B_x_tri_ = L_x_n_*scalar;

	//! A_x_tri_
    construct_L_x_nPlus(t2, index_y);
	A_x_tri_ = L_x_nPlus_*scalar;
	A_x_tri_.addIdMatrix();


	//! F,G
	initlialize_F_G_x(t1, t2, index_y, U_0);

	//! H
	copy_matrixColumnToVector_withoutFirstAndLastElment(Y_0,H_x_t_n_,index_y); // don't care about boudnary condition 
	F_x_t_n_ += H_x_t_n_;
}

void Scheme_CS::calculate_one_step_1 (int index_t_n, int index_y, Matrix& U_0) // X direction:  (1-theta*delta_t*L_x)Y_1      = Y_0 - (theta*deltat_t*L_x)U_0 
{
	double t1      = discret_.get_discret_t(index_t_n);
	double t2      = discret_.get_discret_t(index_t_n+1);
          
	////! Don't care about the boundary condition ! 
	////! Don't care about the boundary condition !
	//if(index_y==0)  //! on the cadre:  v==v_min  --> set directly the value! 
	//{
	//	//cout << "The first boundary condition must can be set as ExtremValue ! " << endl;
	//    for(int index_itr_x=0; index_itr_x<discret_.get_sizeDiscret_x(); ++index_itr_x)
	//	{
	//		Y_1[index_itr_x][index_y] = bc_.bc_Y_L->initialize_cadre(t1, t2,
	//					  											 index_itr_x, index_y,//x, y, 
	//																 U_0);//, U_2);
	//	}
	//}
	//else if(index_y==U_0.cols()-1)  //! on the cadre:  v== v_max --> set directly the value! 
	//{
	//	for(int index_itr_x=0; index_itr_x<discret_.get_sizeDiscret_x(); ++index_itr_x)
	//	{
	//		Y_1[index_itr_x][index_y] = bc_.bc_Y_R->finalize_cadre(t1,t2,
	//															   index_itr_x, index_y,//x, y,
	//															   U_0);//, U_2);
	//	}
	//}
	//else  // solve 1D PDE 
	{
		//! U^n_x: copy index_y'eme column of U --> U_x
		copy_matrixColumnToVector_withoutFirstAndLastElment(U_0,U_1_x,index_y); // without 1st & last element

		//! Construct A,B,H
		initializeScheme_1(t1, t2, index_y, U_0);//, U_1, U_2, scheme_x_itr);  // initialize matrix: A,B,H, need to do : F,G

		//! solve linear Equation:  A*U^{n+1/2} + G^{n+1/2} = B*U^{n}+ F^{n} --> U^{n+1/2} other than the boundary value
		Linear_Equation_Solver::resolve(get_A_x_tri_ref(), // A_x_ will be changed! 
									    get_B_x_tri_ref(),
									    F_x_t_n_,
									    G_x_t_nPlus_,
									    U_1_x,		    // input
									    U_2_x,
										U_x_temp);		    // output

		copy_vectorToMatrixColumn_withoutFirstAndLastElment(Y_1, U_2_x, index_y); //! copy Uplus_x --> index_y'eme column of Uplus

		////! Don't care about the boundary condition ! 
		////! Get the boudnary value of U^{n+1/2}
		//int x_left_border_index  = 0;
		//int x_right_border_index = discret_.get_sizeDiscret_x()-1;

		////! Don't care about the boudnary condition ! 
		//////! if it is the Dirichlet BC or Neumann BC, should use U_2, 
		//////! PDE BC not sure ...  
		//////! Maybe it should be treated before solving the Linear_equation ??? TODO ... 
		//Y_1[0][index_y]            = bc_.bc_X_L->get_bc_extremValue( t1, t2, x_left_border_index,  index_y, U_0);
		//Y_1[U_0.rows()-1][index_y] = bc_.bc_X_R->get_bc_extremValue( t1, t2, x_right_border_index, index_y, U_0);
	}
}

//! scheme_2 direction y:
void Scheme_CS :: initializeScheme_2(double t1, double t2, int index_x, Matrix& U_0)
{
	double scalar = -theta_* (t2-t1); //*delta_t_; //YYYY _changedelta_t_;

	//! B_x_tri_
	construct_L_y_n(t1, index_x);
	B_y_tri_ = L_y_n_*scalar;

	//! A_x_tri_
    construct_L_y_nPlus(t2, index_x);
	A_y_tri_ = L_y_nPlus_*scalar;
	A_y_tri_.addIdMatrix();

	//! F,G
	initlialize_F_G_y(t1, t2, index_x, U_0);

	//! H
	copy_matrixRowToVector_withoutFirstAndLastElment(Y_1,H_y_t_n_,index_x);
	F_y_t_n_ += H_y_t_n_;
}

void Scheme_CS::calculate_one_step_2 (int index_t_n, int index_x, Matrix& U_0)  // Y direction:  (1-theta*delta_t*L_y)U_result = Y_1 - (theta*deltat_t*L_y)U_0
{
	double t1      = discret_.get_discret_t(index_t_n);
	double t2      = discret_.get_discret_t(index_t_n+1);

	//! Boundary condition set to use that obtained by U_0!
	if(index_x==0)   // x==x_min boundary condition
	{
		for(int index_itr_y=0; index_itr_y<discret_.get_sizeDiscret_y(); ++index_itr_y)
		{
			Y_2[index_x][index_itr_y] = bc_->bc_X_L_->initialize_cadre(t1,t2,
																	 index_x,index_itr_y, 
																	 U_0);
		}
	}
	else if(index_x==U_0.rows()-1)  // x=x_max boundary condition
	{
		for(int index_itr_y=0; index_itr_y<discret_.get_sizeDiscret_y(); ++index_itr_y)
		{
			Y_2[index_x][index_itr_y] = bc_->bc_X_R_->finalize_cadre(t1,t2,
																	   index_x,index_itr_y, //x, y,
																	   U_0);//, U_2);
		}
	}
	else  // PDE
	{
		//! UplusHalf_y: copy index_x'eme row of UplusHalf --> UplusHalf_y  
		copy_matrixRowToVector_withoutFirstAndLastElment(U_0, U_1_y, index_x);

		//! Construct A,B,F,G,H
		initializeScheme_2(t1, t2, index_x, U_0);//, U_1, U_2, scheme_y_itr);	

		//! solve linear Equation:  A*U^{n+1} + G^{n+1} = B*U^{n+1/2}+ F^{n+1/2} --> U^{n+1} other than the boundary value
		Linear_Equation_Solver::resolve(  get_A_y_tri_ref(), // A_x_ will be changed! 
										  get_B_y_tri_ref(),
										  F_y_t_n_,
										  G_y_t_nPlus_,
										  U_1_y,           // input
										  U_2_y,
										  U_y_temp);		   // output

		copy_vectorToMatrixRow_withoutFirstAndLastElment(Y_2, U_2_y, index_x); //! copy U_x --> index_x'eme row of U

		int y_left_border_index  = 0;
		int y_right_border_index = discret_.get_sizeDiscret_y()-1;

		Y_2[index_x][0]			   = bc_->bc_Y_L_->get_bc_extremValue(t1, t2, index_x, y_left_border_index,  U_0);//, U_2);
		Y_2[index_x][U_0.cols()-1] = bc_->bc_Y_R_->get_bc_extremValue(t1, t2, index_x, y_right_border_index, U_0);//, U_2);
	}
}







//! ----------------------------------------------------
//! 
//!                   CS part
//! 
//! ----------------------------------------------------
//! scheme_0_prime: Y_0_prime ~ Y_0
void Scheme_CS:: calculate_one_step_0_prime (int index_t_n, Matrix& U_0) 
// Total explicit: Y_0_prime = Y_0 + lambda*L_0(Y_2 - U_0), do it on direction x
{
	//! need boundary condition to calculat the CrossDerivative_
	double t1  = discret_.get_discret_t(index_t_n);
	double t2  = discret_.get_discret_t(index_t_n+1);
	Y_2 -= U_0; 
	caluclate_explicit_cross_derivative(t1, Y_2);
	double scalar = lambda_*(t2-t1); //*delta_t_;
	CrossDerivative_*= scalar;

	//! Don't care about the boundary condition ! 
	add_smallMatrix_to_bigMatrix_withOutBoundary(CrossDerivative_, Y_0);
}


//! Final 
void Scheme_CS :: calculate_one_step(int index_t, Matrix& U_0) //, Matrix& U_1, Matrix& U_2)
{
		double t1	       = discret_.get_discret_t(index_t);         // in fact "t" = tau = T - "t"
		double t2	       = discret_.get_discret_t(index_t+1); 
		//double t_nPlus     = discret_.get_discret_t(index_t+1);

		caluclate_explicit_cross_derivative(t1, U_0);

        //! ------------------------------------------
		//!             Douglas part 
		//! ------------------------------------------
		//! scheme step 0: explicit 
		calculate_one_step_0 (index_t, U_0);

		//! scheme step 1: X direction 
		for(int index_y = 0; index_y<matrix_size.second+2; ++index_y)
		{
			calculate_one_step_1(index_t,  index_y, U_0);
		}

		//! scheme step 2: Y direction 
		for(int index_x=0; index_x<matrix_size.first+2; ++index_x)
		{
			//! ATTENTION Uplus <--> U exchang role
			calculate_one_step_2(index_t,  index_x, U_0);
		}

		//! ------------------------------------------
		//!             CS part 
		//! ------------------------------------------
		//! scheme step 0: explicit 
		calculate_one_step_0_prime (index_t, U_0);

		//! scheme step 1: X direction 
		for(int index_y = 0; index_y<matrix_size.second+2; ++index_y)
		{
			calculate_one_step_1(index_t, index_y, U_0);
		}

		//! scheme step 2: Y direction 
		for(int index_x=0; index_x<matrix_size.first+2; ++index_x)
		{
			//! ATTENTION Uplus <--> U exchang role
			calculate_one_step_2(index_t,  index_x, U_0);
		}

		U_0 = Y_2;

		//Y_0.set_element_to_value(0.0);
		//Y_1.set_element_to_value(0.0);
		//Y_2.set_element_to_value(0.0);
		//U_result.set_element_to_value(0.0);
}