#pragma once
#include <utility>
#include "Scheme_ADI.h"
#include <boost/shared_ptr.hpp>

class Scheme_Yanenko : public Scheme_ADI
{
public:
	//Matrix& U_0;
	mutable Matrix Y;
	mutable Matrix U_result;

	//! constructor & destructor
	Scheme_Yanenko(
					const Discretization& discret,
					PDE_2D_Model& model_discret,
					BoundaryCondition_D2_CONSTPTR bc
				  );

	virtual ~Scheme_Yanenko(){};

	//! Direction x
	void construct_B_x_tri   (double t1, double t2, int index_y);
	void construct_A_x_tri   (double t1, double t2, int index_y);

    //! Direction y
	void construct_B_y_tri    (double t1, double t2, int index_x);
	void construct_A_y_tri    (double t1, double t2, int index_x);

	void initializeScheme_x   (double t1, double t2, int index_y, Matrix& U_0);//, Matrix& U_1, Matrix& U_2, int scheme_x_itr=0);
	void initializeScheme_y   (double t1, double t2, int index_x, Matrix& U_0);//, Matrix& U_1, Matrix& U_2, int scheme_y_itr=0);

	//! one step solver
	void calculate_one_step_x (int index_t_n, int index_y, Matrix& U_0);
	void calculate_one_step_y (int index_t_n, int index_x, Matrix& U_0);

	void  calculate_one_step(int index_t, Matrix& U);

};


typedef boost::shared_ptr<Scheme_Yanenko> Scheme_Yanenko_PTR;
typedef boost::shared_ptr<const Scheme_Yanenko> Scheme_Yanenko_CONSTPTR;