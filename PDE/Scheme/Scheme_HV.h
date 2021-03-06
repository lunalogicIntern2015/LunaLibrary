#pragma once
#include <utility>
#include <PDE/Scheme/Scheme_ADI.h>
#include <boost/shared_ptr.hpp>

class Scheme_HV : public Scheme_ADI
{
public:
	int scheme_itr_size_;
	double theta_;
    double lambda_;

	mutable Matrix Y_0;
	mutable	Matrix Y_1;
	mutable Matrix Y_2;
	mutable Matrix Y_prime; // intermediate variable...
	
public:
	//! constructor & destructor
	Scheme_HV(const Discretization& discret,
			  PDE_2D_Model& model_discret,
			  BoundaryCondition_D2_CONSTPTR bc,
			  double theta );

	virtual ~Scheme_HV(){};


	//! Same as Douglas 
	void initializeScheme_0_x(double t1, double t2, int index_y, Matrix& U_0);
	void initializeScheme_0_y(double t1, double t2, int index_y, Matrix& U_0);
	void initializeScheme_1(double t1, double t2, int index_y, Matrix& U_0);
	void initializeScheme_2(double t1, double t2, int index_y, Matrix& U_0);

	void calculate_one_step_0 (int index_t_n, Matrix& U_0);              // Total explicit: Y_0 = (I+ L_0 + L_1 + L_2)U_0
	void calculate_one_step_1 (int index_t_n, int index_y, Matrix& U_0); // X direction:  (1-theta*delta_t*L_x)Y_1      = Y_0 - (theta*deltat_t*L_x)U_0 
	void calculate_one_step_2 (int index_t_n, int index_x, Matrix& U_0); // Y direction:  (1-theta*delta_t*L_y)U_result = Y_1 - (theta*deltat_t*L_y)U_0

	//! New CS 
	void initializeScheme_0_prime_x(double t1, double t2, int index_y, Matrix& U_0);
	void initializeScheme_0_prime_y (double t1, double t2, int index_y, Matrix& U_0);
	void calculate_one_step_0_prime (int index_t_n, Matrix& U_0);        // Total explicit : Y_0_prime = Y_0 + lambda*L_0(Y_2-U_0)
	void calculate_one_step_1_prime (int index_t_n, int index_y, Matrix& U_0); // X direction:  (1-theta*delta_t*L_x)Y_1      = Y_0 - (theta*deltat_t*L_x)U_0 
	void calculate_one_step_2_prime (int index_t_n, int index_x, Matrix& U_0); // Y direction:  (1-theta*delta_t*L_y)U_result = Y_1 - (theta*deltat_t*L_y)U_0

	void  calculate_one_step(int index_t, Matrix& U);
};

typedef boost::shared_ptr<Scheme_HV> Scheme_HV_PTR;
typedef boost::shared_ptr<const Scheme_HV> Scheme_HV_CONSTPTR;