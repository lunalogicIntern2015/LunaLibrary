#pragma once
#include <utility>
#include <PDE/Matrix/TridiagonalMatrix.h>
#include <PDE/Matrix/Matrix.h>
#include <PDE/Discretization/Discretization.h>
#include <PDE/BoundaryCondition/BoundaryCondition_D2.h>
#include <PDE/PDE_2D_Model.h>
#include <boost/shared_ptr.hpp>
//! --------------------------------------------------------------
//! 
//!     ATTENTION: We suppose time discretization is uniform ! 
//! 
//! --------------------------------------------------------------

class Scheme_ADI
{
public:

	std::pair<int,int> matrix_size;  // x_tilde_size, y_tilde_size

    //! direction x: fist equation
	TridiagonalMatrix L_x_n_;
	TridiagonalMatrix B_x_tri_;
	TridiagonalMatrix L_x_nPlus_;
	TridiagonalMatrix A_x_tri_;

	Matrix F_x_t_n_;          // F = F(boundary condition of B) + H
	Matrix G_x_t_nPlus_;  
	Matrix H_x_t_n_;    // cross-derivative in 1st equation

	//! direction y: second equation
	TridiagonalMatrix L_y_n_;
	TridiagonalMatrix B_y_tri_;
	TridiagonalMatrix L_y_nPlus_;
	TridiagonalMatrix A_y_tri_;

	Matrix F_y_t_n_;         
	Matrix G_y_t_nPlus_;  
	Matrix H_y_t_n_;    // cross-derivative in 2nd equation


	//! always of time t_n ??? (a cause de arbitary boundary condtion of t_{n+1/2} ???)
	Matrix CrossDerivative_;  // of size: (matrix_size.first, matrix_size.second) 

	//SchemasType::SchemasType	   schemasType_; 
	Discretization                 discret_;

	//! keep reference because of polymorphism !
	const PDE_2D_Model&            model_discret_;
	BoundaryCondition_D2_CONSTPTR  bc_;

	//! work place
	mutable Matrix U_1_x;
	mutable Matrix U_2_x;
	mutable Matrix U_x_temp; // for intermediate calculation

	mutable Matrix U_1_y;
	mutable Matrix U_2_y;
	mutable Matrix U_y_temp;  // for intermediate calculation

public:
	//! constructor & destructor
	Scheme_ADI(const Discretization& discret,
				PDE_2D_Model& model_discret,
				BoundaryCondition_D2_CONSTPTR bc
			 );

	virtual ~Scheme_ADI(){};

	//! ****************************************************************
	//!		   STEP_0: Calculate the elemente of Matrix: L_x, L_y
	//! ****************************************************************
	double a(double h, double h_moins) const;
	double b(double h, double h_moins) const;
	double c(double h, double h_moins) const;
	double e(double h, double h_moins) const;
	double f(double h, double h_moins) const;
	double g(double h, double h_moins) const;

	//! ---- ---- Derivative along x ---- ----
	//! ---- Finite Difference
	//! D_1 
	double a_x(double t, int x_index) const;
	double b_x(double t, int x_index) const;
	double c_x(double t, int x_index) const;
	//! D_2
	double e_x(double t, int x_index) const;
	double f_x(double t, int x_index) const;
	double g_x(double t, int x_index) const;
    //! ---- L_s
	double l0_x(double t, int x_index, int y_index) const;
	double l1_x(double t, int x_index, int y_index) const;
	double l2_x(double t, int x_index, int y_index) const;


	//! ---- ---- Derivative along y ---- ----
	//! ---- Finite Difference
	//! D_1 
	double a_y(double t, int y_index) const;
	double b_y(double t, int y_index) const;
	double c_y(double t, int y_index) const;
	//! D_2
	double e_y(double t, int y_index) const;
	double f_y(double t, int y_index) const;
	double g_y(double t, int y_index) const;
    //! ---- L_s
	double l0_y(double t, int x_index, int y_index) const;
	double l1_y(double t, int x_index, int y_index) const;
	double l2_y(double t, int x_index, int y_index) const;



	//! ****************************************************************
	//!		    STEP_1: Explicite Cross Derivative
	//! ****************************************************************
    void  caluclate_explicit_cross_derivative(double t, Matrix& U_t_n) ;
	//! U_n's size              = (matrix_size.first+2, matrix_size.second+2) 
	//! CrossDerivative_'s size = (matrix_size.first, matrix_size.second)

    void get_cross_derivative_x(Matrix& cross_derivative_x, unsigned int index_y)  ;
    void get_cross_derivative_y(Matrix& cross_derivative_y, unsigned int index_x)  ;


	//! ****************************************************************
	//!			 STEP_2: Direction x (y fix)
	//! ****************************************************************
	//! Direction x 
	void construct_L_x_n    (double t1, int index_y)  ; // t_index = t_n_index
	void construct_L_x_nPlus(double t2, int index_y)  ; // t_index = t_nPlus_index

	//! Matrix B_x_trie_ A_x_tri_ should be constructed before calling the following code ! 
	void initlialize_F_G_x    (double t1, double t2, int index_y, Matrix& U_0)  ;
	
	TridiagonalMatrix& get_A_x_tri_ref() {return A_x_tri_;}
	TridiagonalMatrix& get_B_x_tri_ref() {return B_x_tri_;}
	Matrix& get_F_x_t_n_ref()		     {return F_x_t_n_;}
	Matrix& get_G_x_t_nPlus_ref()        {return G_x_t_nPlus_;}


    //! Direction y
	void construct_L_y_n    (double t1, int index_x)  ; // t_index = t_n_index
    void construct_L_y_nPlus(double t2, int index_x)  ; // t_index = t_nPlus_index
	
	void initlialize_F_G_y    (double t1, double t2, int index_x, Matrix& U_0) ;

	const TridiagonalMatrix& get_A_y_tri_ref() const {return A_y_tri_;}
	const TridiagonalMatrix& get_B_y_tri_ref() const {return B_y_tri_;}
	const Matrix& get_F_y_t_n_ref()			 const {return F_y_t_n_;}
	const Matrix& get_G_y_t_nPlus_ref()		 const {return G_y_t_nPlus_;}

	//! calculate one step
	virtual void calculate_one_step(int index_t, Matrix& U_0)  =0;
};

typedef boost::shared_ptr<Scheme_ADI> Scheme_ADI_PTR;
typedef boost::shared_ptr<const Scheme_ADI> Scheme_ADI_CONSTPTR;