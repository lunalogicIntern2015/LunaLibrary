#include <iostream>
#include <sstream>
#include <fstream>

#include "PDE_ADI_Solver.h"
#include "HestonModel.h"
#include "useful_function.h"
#include "Params.h"

using namespace std;

//! ---- OUTPUT ---- 
const bool DEBUG_printing_result       = false;
const bool DEBUG_printing_result_final = false;	

//const string DEBUG_output_path = "C:\\Documents and Settings\\n87206\\My Documents\\Visual Studio 2005\\Projects\\PDE\\PDE_D2\\Output\\";

PDE_ADI_Solver::PDE_ADI_Solver(Scheme_ADI& scheme)
		      : 
              scheme_(scheme),
			  U      (scheme_.matrix_size.first+2,scheme_.matrix_size.second+2,0)
{ 
     cout << "matrix U's size = (" <<  scheme_.matrix_size.first << "," << scheme_.matrix_size.second << ")"  << endl;
}

 
PDE_ADI_Solver :: ~PDE_ADI_Solver(){} 


void PDE_ADI_Solver::Initialize_U()
{
	for(int i=0; i< scheme_.matrix_size.first+2; i++)
	{
		for(int j=0; j<scheme_.matrix_size.second+2; j++)
		{
			U[i][j] = scheme_.bc_->bc_I_->get_inital_bc(scheme_.discret_.get_discret_x(i),scheme_.discret_.get_discret_y(j));  //! v is not useful
		}
	}
}


//! direct calculation
void PDE_ADI_Solver::solve_PDE()
{
	Discretization& discret   = scheme_.discret_;
	BoundaryCondition_D2_CONSTPTR bc  = scheme_.bc_;

	//! t_0
	Initialize_U(); 

	//! ---- OUTPUT 1 ----
	if(DEBUG_printing_result == true)
	{
		std::string  file_name_U_original = DEBUG_output_path + "Original.csv";
		print_result(file_name_U_original, discret.get_discret_x(),discret.get_discret_y(),  U);
	}

	//! evaluate U_i --> U_{i+1}
    for(int index_t =0; index_t <discret.get_sizeDiscret_t()-1; ++index_t)  // time ! 
	{
		double t_n	       = discret.get_discret_t(index_t);         // in fact "t" = tau = T - "t"
		double t_nPlus     = discret.get_discret_t(index_t+1);

		cout.precision(6);
		
		//cout << "PDE solving t: [" << t_n << "," << t_nPlus << "], ~~~~ " << index_t/(double)(discret.get_sizeDiscret_t()-1)<< "%"<< endl;
		
		if(scheme_.model_discret_.get_model_discret_type() != "SLV")  // normal case!
		{
			scheme_.calculate_one_step(index_t, U);
		}
		else
		{
			int num_iteration = 2;
			for(int itr = 0; itr<num_iteration; ++itr)
			{
				scheme_.model_discret_.calculate_L(index_t, U);
				scheme_.calculate_one_step(index_t, U);
			}
		}
		//! if SLV need to recursion! 
	}
}

