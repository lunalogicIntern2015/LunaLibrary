//! Question: when Feller condition is satisfied, and S0=100, v0=0.04, why the sum proba is not 1? (only about 0.88)
//! when S0=100, v0=0.2, much better (sum proba = 0.97), but cannot improved to 1, which is different from 1D case ... ??????

//! For changing the code: check
//! PDE:        18.7140820773669
//! Analytical: 18.7998851019194



//! ---- ---- ---- ---- ---- ---- ---- ---- 
//!				Forward PDE 
//! ---- ---- ---- ---- ---- ---- ---- ---- 
//! Rewrite class discretization: about the shift part, too bad ...
//! 
//! 
//! Observation 1.  
//! Fwd PDE: non uniform discretizaton on t is better than uniform: they have the same sum-proba, but 
//!          when we have a closer look, much less negative proba,
//!          One Example: (nonuniform_scale_param = 1.0/50)
//!          mixed non-uniform+uniform: sump_proba = 0.968921, sump_proba_negative = -0.0340168,  sump_proba_positive = 1.00294 (scale_param = 1/50)
//!          uniform case : sump_proba = 0.972173,  sump_proba_negative = -0.0514,    sump_proba_positive = 1.02359   total_error = 0.075
//!			 non uniform  : sump_proba = 0.968848,  sump_proba_negative = -0.0000195, sump_proba_positive = 0.968868, total_error = 0.031   (scale_param = 1/50)
//! 		 non uniform  : sump_proba = 0.9689537, sump_proba_negative = -0.05139,   sump_proba_positive = 1.02093,  total_error = 0.0723   (scale_param = 1/5)
//! 		 non uniform  : sump_proba = 0.971447,  sump_proba_negative = -0.056433,  sump_proba_positive = 1.02788,  total_error = 0.084   (scale_param = 1)
//! 		 non uniform  : sump_proba = 0.968853,  sump_proba_negative = -0.0000028,  sump_proba_positive = 0.968855,  total_error = 0.0311474   (scale_param = 1)


//! Observation 2.
//! une observaton tres bizzar, when s,v grid extreamly small, and time step is big, then: theta = 0.5 sheme seems better than theta = 1 schem,
//! in the sense of the negative probability.
//!		in fact when time step == 1, theta=1 is better than theta == 0.5, but once the time step is more than 2 or 3, theta =0.5 can achieve 0
//!     negative proba, but theta==1 shceme has always negative proba ... (even in the case crossing item = 0) Theta=0.5 seems can automatically
//!     correct the negative probability.
//! Don't understand why ... ... 


//! Observation 3.
//! When use theta=0.5 scheme, can achieve negative proba =0 when crossing item =0, but when crossing item != 0 there will 
//! always be some negative probability.


//! Observation 4 (bizzar) ??? ??? ??? ??? ??? 
//! when theta = 0.5 seems to be more stable than theta  =1; 
//! Parameters used: 
//! t discretization/scale:  30,  1/50000
//! x discretization/scale:  100, K/5
//! y discretization/scale:  100, 1
//const double Heston_params_kappa = 1.5;
//const double Heston_params_theta = 0.04;
//const double Heston_params_sigma = 0.3;
//const double Heston_params_rho   = -0.5;
//const double Heston_params_r     = 0.025;
//const double Option_params_K = 100.0;
//const double Option_params_T = 1.0; 
//const double fwd_PDE_S0 = 100.0;
//const double fwd_PDE_v0 = 0.2;


//! observation 5
//! Heston_param_theta = 0.04
//! When v_0 = 0.04 or 0.2, need different discretization parameters to achieve 0.99 sum_positive_proba, and 10E-4 sum_negative_proba
//! v_0 = 0.2:  mixed discretization on y: 
//!     nonuniform --> t: 50,1/1000 
//!     mixed on y --> nonuniform: 50, 1/10000 (or 1/1000)   uniform-->y:50
//!     <>
//! v_0 = 0.04: only nonuniform on y (if use mixed discretization, the sum_positive_proba will decrease .. WHY ???): 
//!     nonuniform --> t: 25, 1/1000
//!     nonuniform --> y: 50, 1/10000  
//!    <don't shift>


//! observation 6, 
//! Another pricing formula, need to *= (2/3) ??? 




//! Fwd PDE

#include <iostream>
#include <fstream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>

#include "Discretization.h"
#include "BoundaryCondition_Heston_Put.h"
#include "BoundaryCondition_D2_density.h"

#include "PDE_ADI_Solver.h"
//#include "PDE_ADI_Extrapolation_Solver.h"

#include "Scheme_ADI.h"
#include "Scheme_Yanenko.h"  // 1
#include "Scheme_Douglas.h"  // 2 
#include "Scheme_CS.h"       // 3
#include "Scheme_MCS.h"      // 4
#include "Scheme_HV.h"       // 5

//#include "HestonModel.h"
#include <PDE/PDE_Model_Heston.h>
#include "Model.h"
#include "Heston_analytical_pricer.h"
#include "useful_function.h"

#include <PDE/Matrix/Matrix.h>
#include "Params.h"
#include <fstream>

#include "Interpolator.h"
/*
Heston Model: 
  \frac{dS}{S} = r*dt + sqrt(v)*dW^1
  dv = \kappa(\theta - v)*dt + sigma*sqrt(v)*dW^2
  <dW^1,dW^2> = \rho*dt
*/

//! some constat about the stability: 
//! 1. When the sinh non-unifor scheme is used, n_t_discret = 100; n_x_discre = 200, n_y_discret = 50:
//! Yanenko is stable, Douglas (1) is stable,  CS (1, 1/2) is stable
//! Douglas(1/2) not stable !  CS (1/2, 1/2) not stable ! 


const int    scheme_type           = 3;   // 1: Yanenko, 2:Douglas, 3:CS, 4:MCS, 5:HV
const double scheme_douglas_theta  = 0.5; // used for scheme: 2,3,4,5
const double scheme_douglas_lambda = 0.5; // used for scheme: 3,4



const bool  DEBUG_test_PDE_solver_extrapolation   = false; // Don't work ... ! What to do for 2D extrapolation ???

const bool  DEBUG_PDE_test_only_PDE				  = true; 
const bool  DEBUG_PDE_test_PDE_vs_analytical	  = !DEBUG_PDE_test_only_PDE;

const bool  DEBUG_PDE_error_print_to_file		  = true;  //! for all grid (S,v)
const bool  DEBUG_PDE_error_print_to_screen		  = false;  //! for only one grid :)

const double CONFIG_S_max = Option_params_K*5;//Option_params_K*8.0;
const double CONFIG_v_max = 0.4;//5.0;

const int    CONFIG_multiplicator = 1;
const int    CONFIG_sizeDiscretization_t = 100*CONFIG_multiplicator;
const int    CONFIG_sizeDiscretization_x = 100*CONFIG_multiplicator; //Attention: when compare the precision of uniform and non-uniform result, we need a grid at S = 100 !
//const int    CONFIG_sizeDiscretization_y = 50*CONFIG_multiplicator; 
const bool   OPTIMIZATION_if_nonuniform_sinh_grid_t = true;   
const bool   OPTIMIZATION_if_shfit_grid_t           = false;   
const bool   OPTIMIZATION_if_nonuniform_sinh_grid_x = true;   // nonuniform discretization on space: x-direction
const bool   OPTIMIZATION_if_shfit_grid_x           = true;   // OPTIMIZATION_if_nonuniform_sinh_grid_x;
const bool   OPTIMIZATION_if_nonuniform_sinh_grid_y = true;   // 0onuniform discretization on space: y-direction
const bool   OPTIMIZATION_if_shfit_grid_y           = false;
const double OPTIMIZATION_nonuniform_t_sinh_scale_param = 1/10000.0;
//const double OPTIMIZATION_nonuniform_x_sinh_scale_param = Option_params_K/5; // Heston pricing
const double OPTIMIZATION_nonuniform_x_sinh_scale_param = 1/1.0; // Fwd PDE
const double OPTIMIZATION_nonuniform_y_sinh_scale_param = 1/10000.0;

//! ------------------------------------------------------------`
//! In current setting: S_max = 800, v_max = 5, n_t = 100, n_x = 160, n_y=100
//! max error for interval: [0,1.5K] x [0,1]
//! (u,u) = -0.66  at position (100,0)
//! (u,n) = -0.12  at position (100,0.95)
//! (n,n) = -0.11  at position (100,0.95)
//!
//!  max error for interval: [0,1.5K] x [0,0.1]
//! (u,u) = -0.66  at position (100,0)
//! (u,n) = -0.12  at position (100,0.00069)
//! (n,n) = -0.026 at position (100,0.097)
//! ------------------------------------------------------------

////! BlackModel: dS/S = r*dt + sigma*dW  (in Heston sigma ~ theta)
//double estimat_x_max(double S0, double v_max, double T)  // proba of S_T between [0,x_max]
//{
//	double r    = Heston_params_r;
//	double vol  = v_max;
//
//	cout << r << "    ----    " << vol << endl;
//	
//	//G<4:    0.99996
//	//G<4.5 : 0.999996
//	//G<5 :   0.9999997
//	double a  = 4.5;
//	return S0*exp( a*vol*sqrt(T) + (r-vol*vol/2)*T );
//}


Interpolator_D2 test_Heston_PDE_solution(int CONFIG_sizeDiscretization_y)
{	
	std::cout << "Model and Product Direction_Parameters" << std::endl;

	//! Heston Direction_Parameters
	double kappa = Heston_params_kappa;
	double theta = Heston_params_theta;
	double sigma = Heston_params_sigma;
	double rho   = Heston_params_rho;
	double r     = Heston_params_r;

	//! Product Direction_Parameters
	double K = Option_params_K;
	double T = Option_params_T; 

	//! FD param
	double S_max = CONFIG_S_max;  // S_max = 2*K is not efficient! 
	double v_max = CONFIG_v_max;

	//! resume: the biggest error is the time! 
	int sizeDiscretization_t = CONFIG_sizeDiscretization_t + 1;  
	int sizeDiscretization_x = CONFIG_sizeDiscretization_x + 1; 
	int sizeDiscretization_y = CONFIG_sizeDiscretization_y + 1;

 //   //! ma model
	//MaModel ma_model(Ma_params_alpha,Ma_params_sigma,Ma_params_rho);
	//Ma_Model_discret ma_model_discret(ma_model);

	//! Range and discretization
	Range range_t(0.0,  T);        
	Range range_x(0.0,  S_max);   
	Range range_y(0.0,  v_max);  


	//! non-uniform 
		//! direction t
	int t_discretization_size = 1;
	std::vector<bool>   if_nonuniform_t(t_discretization_size);//              = OPTIMIZATION_if_nonuniform_sinh_grid_x;
    std::vector<bool>   if_shift_t(t_discretization_size); //                  = OPTIMIZATION_if_shfit_grid_x;
	std::vector<double> t_nonuniform_center_c(t_discretization_size);//        = K; 
	std::vector<double> t_nonuniform_scale_param_c(t_discretization_size);//   = K/5.0; 
	std::vector<double> t_shifting_center(t_discretization_size); //            = K;

	if_nonuniform_t[0]              = OPTIMIZATION_if_nonuniform_sinh_grid_t;
    if_shift_t[0]                   = OPTIMIZATION_if_shfit_grid_t;
	t_nonuniform_center_c[0]        = 0.0; 
	t_nonuniform_scale_param_c[0]   = OPTIMIZATION_nonuniform_t_sinh_scale_param; //1.0/500.0; 
	t_shifting_center[0]            = -999999.9;

	if(if_fwd_PDE == true)
	{
		////! add second time
  //     	if_nonuniform_t.push_back(false);             //           = true;//OPTIMIZATION_if_nonuniform_sinh_grid_t;
		//if_shift_t.push_back(false)    ;              // = OPTIMIZATION_if_shfit_grid_t;
		//t_nonuniform_center_c.push_back(0.0);		  //        = 0.0; 
		//t_nonuniform_scale_param_c.push_back(1.0/50); //   = 1.0/25.0; 
		//t_shifting_center.push_back(-99999.9);        //            = -999999.9;
	}

		//! direction x
	int x_discretization_size = 1;
	std::vector<bool>   if_nonuniform_x(x_discretization_size);//              = OPTIMIZATION_if_nonuniform_sinh_grid_x;
    std::vector<bool>   if_shift_x(x_discretization_size); //                  = OPTIMIZATION_if_shfit_grid_x;
	std::vector<double> x_nonuniform_center_c(x_discretization_size);//        = K; 
	std::vector<double> x_nonuniform_scale_param_c(x_discretization_size);//   = K/5.0; 
	std::vector<double> x_shifting_center(x_discretization_size); //            = K;

	if_nonuniform_x[0]              = OPTIMIZATION_if_nonuniform_sinh_grid_x;
    if_shift_x[0]                   = OPTIMIZATION_if_shfit_grid_x;
	x_nonuniform_center_c[0]        = K; 
	x_nonuniform_scale_param_c[0]   = K/5.0; 
	x_shifting_center[0]            = K;

	if(if_fwd_PDE == true)
	{
		x_nonuniform_center_c[0]      = fwd_PDE_S0; 
		x_nonuniform_scale_param_c[0] = OPTIMIZATION_nonuniform_x_sinh_scale_param; 
		x_shifting_center[0]          = fwd_PDE_S0;
	}

		//! direction y
	int y_discretization_size = 1;
	std::vector<bool>   if_nonuniform_y(y_discretization_size);		    
    std::vector<bool>   if_shift_y(y_discretization_size);				
	std::vector<double> y_nonuniform_center_c(y_discretization_size);   
	std::vector<double> y_nonuniform_scale_param_c(y_discretization_size); 
	std::vector<double> y_shifting_center(y_discretization_size);          

	if_nonuniform_y[0]		      = OPTIMIZATION_if_nonuniform_sinh_grid_y;
    if_shift_y[0]				  = OPTIMIZATION_if_shfit_grid_y; //! for Fwd PDE
	y_nonuniform_center_c[0]      = 0.0;     
	y_nonuniform_scale_param_c[0] = OPTIMIZATION_nonuniform_y_sinh_scale_param;  // V_max/500 is the experienced choice! 
	y_shifting_center[0]          = -99999.9;

	if(if_fwd_PDE == true)
	{
		//! first dimension: 
		if_nonuniform_y[0]			  = OPTIMIZATION_if_nonuniform_sinh_grid_y; // OPTIMIZATION_if_nonuniform_sinh_grid_y;
		if_shift_y[0]				  = OPTIMIZATION_if_shfit_grid_y; // OPTIMIZATION_if_shfit_grid_y; //! for Fwd PDE
		y_nonuniform_center_c[0]      = 0.0;//Heston_params_theta;     
		y_nonuniform_scale_param_c[0] = OPTIMIZATION_nonuniform_y_sinh_scale_param; // OPTIMIZATION_nonuniform_y_sinh_scale_param; //1.0/10000.0;  // V_max/500 is the experienced choice! 
		y_shifting_center[0]          = fwd_PDE_v0;

		////! second dimension:
		//if_nonuniform_y.push_back(false);
		//if_shift_y.push_back(true);					                  
		//y_nonuniform_center_c.push_back(fwd_PDE_v0);			      
		//y_nonuniform_scale_param_c.push_back(1.0/*OPTIMIZATION_nonuniform_y_sinh_scale_param*/);
		//y_shifting_center.push_back(fwd_PDE_v0); 

		////! third
		//if_nonuniform_y.push_back(false);
		//if_shift_y.push_back(true);					                  
		//y_nonuniform_center_c.push_back(Heston_params_theta);			      
		//y_nonuniform_scale_param_c.push_back(1.0/5000.0);
		//y_shifting_center.push_back(fwd_PDE_v0);
	}

	Discretization discret( range_t, sizeDiscretization_t,  
							range_x, sizeDiscretization_x,  
							range_y, sizeDiscretization_y,
							if_nonuniform_t,  t_nonuniform_center_c,  t_nonuniform_scale_param_c,
							if_shift_t,    t_shifting_center,
							if_nonuniform_x,  x_nonuniform_center_c,  x_nonuniform_scale_param_c,
							if_shift_x,    x_shifting_center,
							if_nonuniform_y,  y_nonuniform_center_c,  y_nonuniform_scale_param_c,
							if_shift_y,    y_shifting_center);  


	//! Heston Model
	HestonModel model(r,sigma,kappa,theta,rho);
	//Heston_Model_discret heston_model_discret(model);
	//! Heston discret Model
	PDE_Model_Heston heston_model_discret(model,discret);


    //! Term Structure, Payoff 	&&  Boundary Condition
	//! Don't want to use this one, it will make error in the code ! 
	//! Because Heston model calls directly r() of Heston model :) 
	Term_Structure* cts  = NULL; //new Const_Term_Structure(-999999);  
	PayoffYuan* payoff   = new Put(K);

	//! Backward PDE
	//BoundaryCondition_Heston_Put bc(cts, payoff, model, discret);
	//! Forward PDE
	BoundaryCondition_D2_density_CONSTPTR bc (new BoundaryCondition_D2_density (model, discret,fwd_PDE_S0,fwd_PDE_v0, cts, payoff ));

	//! --------------------------------------------------
	//! 
    //!           Test PDE solver 
	//! 
	//! ---------------------------------------------------
	//! new Scheme class 
	Scheme_ADI* scheme_new_ptr = NULL;
	if(scheme_type ==1)
	{
		scheme_new_ptr= new Scheme_Yanenko(
			discret,
			heston_model_discret,
			bc
			);
	}
	else if(scheme_type ==2)
	{
		scheme_new_ptr = new Scheme_Douglas(
			discret,
			heston_model_discret,
			bc,
			scheme_douglas_theta
			);
	}
	else if (scheme_type == 3)
	{
		scheme_new_ptr = new Scheme_CS(
			discret,
			heston_model_discret,
			bc,
			scheme_douglas_theta,
			scheme_douglas_lambda
			);
	}
	else if(scheme_type == 4)
	{
		scheme_new_ptr = new Scheme_MCS(
			discret,
			heston_model_discret,
			bc,
			scheme_douglas_theta,
			scheme_douglas_lambda
			);
	}
	else if(scheme_type == 5)
	{
		scheme_new_ptr = new Scheme_HV(
			discret,
			heston_model_discret,
			bc,
			scheme_douglas_theta
			);
	}
	else
	{
	    throw ("Error in Test_PDE.cpp, scheme_type  does not take valid value");
	}

	PDE_ADI_Solver   pde_solver(*scheme_new_ptr);

	//! --------------------------------------------------
	//! 
    //!           Test PDE solver extrapolation
	//! 
	//! ---------------------------------------------------
	//PDE_ADI_Extrapolation_Solver pde_extrapolation_solver( scheme_type,
	//													   discret,
	//													   heston_model_discret,
	//													   // --- For constructing: stupid bc (Bad Bad design !!!)
	//													   cts, 
	//													   payoff,
	//													   model,
	//													   // ----
	//													   scheme_douglas_theta, 
	//													   scheme_douglas_lambda);

	clock_t t1 = clock();
	//if(DEBUG_test_PDE_solver_extrapolation == true)
	//{
	//	pde_extrapolation_solver.solve_PDE();
	//}
	//else
	{
 		pde_solver.solve_PDE();	
	}
	clock_t t2 = clock();
	cout << "time = " << (double)(t2-t1)/CLOCKS_PER_SEC << " seconds"<< endl;

	////! Fwd PDE
	////! find initial dirac point
	//int dirac_x_index = -1;
	//int dirac_y_index = -1;
	//vector<double> discret_x_tilde = discret.get_discret_x_tilde();
	//vector<double> discret_y_tilde = discret.get_discret_y_tilde();
	//for(unsigned int i=0; i<discret_x_tilde.size(); ++i)
	//{
	//	double x = discret_x_tilde[i];
	//	if(fabs(x-fwd_PDE_S0)<0.000000001)
	//	{
	//		dirac_x_index = i;  
	//		break;
	//	}
	//}
	//for(unsigned int j=0; j<discret_y_tilde.size(); ++j)
	//{
	//	double y = discret_y_tilde[j];
 //       if(fabs(y-fwd_PDE_v0)<0.000000001)
	//	{
	//		dirac_y_index = j;
	//		break;
	//	}
	//}

	//
	//double inital_proba = 1.0/4
	//	                  *( discret.get_delta_x(dirac_x_index)+discret.get_delta_x(dirac_x_index+1) )
	//					  *( discret.get_delta_y(dirac_y_index)+discret.get_delta_y(dirac_y_index+1) );
	//cout << "inital_proba  = " << inital_proba  << endl;

	double inital_proba = 1.0;

	//! check sum of probability:
	Matrix Fwd_PDE_result (pde_solver.get_result());
	Fwd_PDE_result /= inital_proba;

	double sum_proba          = 0.0;
	double sum_proba_positive = 0.0;
	double sum_proba_negative = 0.0;
	for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
	{
		int index_i = i+1;
		for( int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
		{
			int index_j = j+1;
            double density = Fwd_PDE_result[index_i][index_j]; // then it is the density
			double coeff   = 1.0/4
		                    *( discret.get_delta_x(i)+discret.get_delta_x(i+1) )
						    *( discret.get_delta_y(j)+discret.get_delta_y(j+1) );

			double proba   = coeff*density;  
			if(density>0.0)
			{
				sum_proba_positive += coeff*density;  
			}
			else
			{
			   sum_proba_negative  += coeff*density;  
			}
		}
	}
	sum_proba = sum_proba_positive + sum_proba_negative;
	cout << " ---- |||| ---- |||| Sum_Proba          = " << sum_proba << endl;  //<< ",  sum_proba_positive = " << sum_proba_positive << " + " << "sum_proba_negative = " << sum_proba_negative << endl;
	cout << " ---- |||| ---- |||| Sum_Proba_positive = " << sum_proba_positive << endl;
	cout << " ---- |||| ---- |||| Sum_Proba_negative = " << sum_proba_negative << endl;
	cout << " ---- |||| ---- |||| Sum_total_error    = " << fabs(sum_proba_positive - 1) + fabs(sum_proba_negative) << endl;


	//! analytical price
	int num_K = 17;
	double step = 0.05;
	//double K_start = Option_params_K*(num_K-1)/

	vector<double> v_K(num_K);
	vector<double> v_analytical_price(num_K);

	cout << "strike vector: "<< endl;
	cout << Option_params_K << endl;
	v_K[0] = Option_params_K;
	for(unsigned int i=1; i<v_K.size()-1; i+=2)
	{
		v_K[i]   =  Option_params_K *(1+i*step);
		v_K[i+1] =  Option_params_K *(1-i*step);
		cout << v_K[i] << endl;
		cout << v_K[i+1] << endl;
	}

	for(unsigned int i=0; i<v_K.size(); ++i)
	{
		Heston_analytical_pricer pricer_ana(v_K[i]);
		v_analytical_price[i] = pricer_ana.heston_put(fwd_PDE_S0,fwd_PDE_v0, 1024*10);

	}

	vector<double> v_numerical_price(num_K);
	//! first numerical integration price
	for(unsigned int kk=0; kk<v_K.size()-1; ++kk)
	{
		    double KK = v_K[kk];
			double sum_price          = 0.0;
			double sum_price_positive = 0.0;
			double sum_price_negative = 0.0;
			for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
			{
				int index_i = i+1;
				double payoff = 1;
				double x = discret.get_discret_x(index_i);
				if(x<KK)
					payoff = KK - x;
				else
					payoff = 0.0;

				for( int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
				{
					int index_j = j+1;
					double density = Fwd_PDE_result[index_i][index_j]; // then it is the density
					double coeff   = 1.0/4
									*( discret.get_delta_x(i)+discret.get_delta_x(i+1) )
									*( discret.get_delta_y(j)+discret.get_delta_y(j+1) );

					double proba   = coeff*density;  
					if(density>0.0)
					{
						sum_price_positive += coeff*density*payoff;  
					}
					else
					{
					   sum_price_negative  += coeff*density*payoff;  
					}
				}
			}
			sum_price_positive *= exp(-Option_params_T*Heston_params_r);
			sum_price_negative *= exp(-Option_params_T*Heston_params_r);
			sum_price = sum_price_positive + sum_price_negative;
			v_numerical_price[kk] = sum_price;
			cout << "1st numerical price = " << sum_price << endl;
			//cout << "price_numerical_positive  = " << sum_price_positive << endl;
			//cout << "price_numerical_negative  = " << sum_price_negative << endl;


			//!second numerical integration price
			vector<double> x_grid(discret.get_discret_x_tilde());
			vector<double> y_grid(discret.get_discret_y_tilde());

			vector<double> x_value(discret.get_sizeDiscret_x_tilde(),0.0);
			vector<double> y_value(discret.get_sizeDiscret_y_tilde(),0.0);

			for(int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
			{
				int index_i = i+1;
				double payoff = 1;
				//double x = discret.get_discret_x(index_i);
				//if(x<KK)
				//	payoff = KK - x;
				//else
				//	payoff = 0.0;

				for(int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
				{
					 int index_j = j+1;
					 y_value[j] = Fwd_PDE_result[index_i][index_j];
				}
				Interpolator interp_y(y_grid,y_value);
				x_value[i] = interp_y.calculate_integral_approximation()*payoff;
			}
			Interpolator interp_x(x_grid,x_value);
			double second_analytical_price = interp_x.calculate_integral_approximation();
			//second_analytical_price *= exp(-Option_params_T*Heston_params_r);
			cout << "2nd numerical price = " << second_analytical_price << endl;
			

			//! third numerical integration price
			double third_sum_price          = 0.0;
			vector<double> v_sum_x(discret.get_sizeDiscret_x_tilde(),0.0);
			vector<double> v_sum_x_positive(discret.get_sizeDiscret_x_tilde(),0.0);
			vector<double> v_sum_x_negative(discret.get_sizeDiscret_x_tilde(),0.0);
			for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
			{
				int index_i = i+1;
				double payoff = 1;
				//double x = discret.get_discret_x(index_i);
				//if(x<KK)
				//	payoff = KK - x;
				//else
				//	payoff = 0.0;
				double sum_x          = 0.0;
				double sum_x_positive = 0.0;
				double sum_x_negative = 0.0;
				for( int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
				{
					int index_j = j+1;
					double density = Fwd_PDE_result[index_i][index_j]; // then it is the density
					double coeff   = ( discret.get_delta_y(j)+discret.get_delta_y(j+1) )/2;

					double proba   = coeff*density;  
					if(density>0.0)
					{
						sum_x_positive += coeff*density*payoff;  
					}
					else
					{
					   sum_x_negative  += coeff*density*payoff;  
					}
					sum_x = sum_x_positive + sum_x_negative;
				}
				v_sum_x_positive[i] = sum_x_positive;
				v_sum_x_negative[i] = sum_x_negative;
				v_sum_x[i]          = sum_x;
			}

			double hehe = 0.0;
			for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
			{
				double interval = ( discret.get_delta_x(i)+discret.get_delta_x(i+1) )/2.0;
				hehe += v_sum_x[i]*interval;
			}
			cout << "3rd numerical price = " << hehe << endl;
	}
	Interpolator_D2 interp_2d(discret.get_discret_x(), discret.get_discret_y(), Fwd_PDE_result);

			////sum_price_positive *= exp(-Option_params_T*Heston_params_r);
			////sum_price_negative *= exp(-Option_params_T*Heston_params_r);
			//sum_price = sum_price_positive + sum_price_negative;
			//cout << "price_numerical  = " << sum_price << endl;
			//cout << "price_numerical_positive  = " << sum_price_positive << endl;
			//cout << "price_numerical_negative  = " << sum_price_negative << endl;
	


	//! Print one price
	if(DEBUG_PDE_error_print_to_screen == true)
	{
		cout.precision(15); 

		int S_index = 25;
		int v_index = 10;

		double S_check = discret.get_discret_x(S_index);
		double v_check = discret.get_discret_y(v_index);

		//Matrix PDE_result       = pde_solver.get_result();
		Matrix PDE_result(1,1,0);
		//if(DEBUG_test_PDE_solver_extrapolation == true)
		//{
		//	PDE_result = pde_extrapolation_solver.get_result();
		//}
		//else
		{
		    PDE_result = pde_solver.get_result();
		}

		Heston_analytical_pricer analytical_pricer;
		vector<double> discret_x;
		discret_x.push_back(S_check);
		
		cout << "PDE's price for So=" << S_check << ", v=" << v_check <<  "   =   " << PDE_result[S_index][v_index] << endl;

		vector<double> discret_y;
		discret_y.push_back(v_check);
		Matrix analytical_price = analytical_pricer.heston_put_v(discret_x,discret_y, 1024*100);
		cout.precision(15); 
		cout << "analytical formula's price = " << analytical_price[0][0]<< endl;
	}

	if (DEBUG_PDE_error_print_to_file == true)
	{
		std::stringstream ss;
		ss << CONFIG_sizeDiscretization_y;
		string s = ss.str();
		std::string PDE_error_output_file_name = DEBUG_output_path + s + "PDE_error.csv";

		//! get PDE result
		//Matrix PDE_result       = pde_solver.get_result();
		//Matrix PDE_result(1,1,0);
		//if(DEBUG_test_PDE_solver_extrapolation == true)
		//{
		//	PDE_result = pde_extrapolation_solver.get_result();
		//}
		//else
		//{
		//    PDE_result = pde_solver.get_result();
		//}

		Matrix PDE_result = Fwd_PDE_result;

		//! get analytical result 
		Heston_analytical_pricer analytical_pricer;

		vector<double> discret_x(discret.get_discret_x());
		vector<double> discret_y(discret.get_discret_y());

		//! Backward PDE
		//double Error_S_min = 0.5 * Option_params_K;
		//double Error_S_max = 1.5 * Option_params_K;
		//double Error_v_min = 0;
		//double Error_v_max = 1.0;

		//! Fwd PDE
		double Error_S_min = 0;
		double Error_S_max = S_max; // 1.5 * Option_params_K;
		double Error_v_min = 0;
		double Error_v_max = v_max; //1.0;


		vector<double> x_vector_value;  //! For analytical price
		vector<int>    x_vector_index;  
		vector<double> y_vector_value;  //! For PDE price (index) 
		vector<int>    y_vector_index;   

		for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
		{
			double x = discret_x[i];
		    if(x>=Error_S_min-0.000000001 && x <= Error_S_max+0.000000001)
			{
				x_vector_value.push_back(x);
				x_vector_index.push_back(i);
			}
		}

		for(unsigned int i=0; i<discret_y.size(); ++i)
		{
			double y = discret_y[i];
		    if(y>=Error_v_min-0.000000001 && y <= Error_v_max+0.000000001)
			{
				y_vector_value.push_back(y);
				y_vector_index.push_back(i);
			}
		}
		
		Matrix analytical_price(1,1,0.0);
		if(DEBUG_PDE_test_PDE_vs_analytical == true)
		{
			analytical_price  = analytical_pricer.heston_put_v(x_vector_value,y_vector_value);
		}
		Matrix PDE_error((int)x_vector_value.size(), (int)y_vector_value.size(),-9999.0);  
		for(unsigned int i=0; i<x_vector_value.size(); ++i)
		{
			for(unsigned int j=0; j<y_vector_value.size(); ++j)
			{
			    if(DEBUG_PDE_test_PDE_vs_analytical == true)
				{
				    PDE_error[i][j] = PDE_result[x_vector_index[i]][y_vector_index[j]] - analytical_price [i][j];
				}
				else if (DEBUG_PDE_test_only_PDE == true)
				{
				    PDE_error[i][j] = PDE_result[i][j];
				}
			}
		}
		print_result(PDE_error_output_file_name , x_vector_value, y_vector_value, PDE_error);

		//! find max error
		if(DEBUG_PDE_test_PDE_vs_analytical == true)
		{
			double max_error_value   = 0.0;
			int    max_error_index_S = 0;
			int    max_error_index_v = 0;
			for(unsigned int i=0; i<x_vector_value.size(); ++i)
			{
				for(unsigned int j=0; j<y_vector_value.size(); ++j)
				{
					if(abs(PDE_error[i][j]) > abs(max_error_value))
					{
						max_error_value = PDE_error[i][j];
						max_error_index_S = i;
						max_error_index_v = j;
					}
				}
			}
			//!!!!! I donot know why but Error !!!!  -- Max error = 17
			cout << "---- ---- ---- ---- Max_error_value is at (" << x_vector_value[max_error_index_S]<< ","<< y_vector_value[max_error_index_v] << ")   = " <<  max_error_value  << endl;
		}
	}

	for(unsigned int i=0; i<v_K.size(); ++i)
	{
		cout << "v_K[" << i << "] = "<< v_K[i] << " 's price: " << v_analytical_price[i] << " ~ " << v_numerical_price[i] << endl;
	}

	delete cts;
	delete payoff;
	//delete scheme_ptr;
	delete scheme_new_ptr;

	return interp_2d;
}


void test_Heston_PDE()
{	
	cout << "Parameters used: " << endl;
	cout << "CONFIG_sizeDiscretization_t = "<< CONFIG_sizeDiscretization_t << endl;
	cout << "CONFIG_sizeDiscretization_x = "<< CONFIG_sizeDiscretization_x << endl;
	//cout << "CONFIG_sizeDiscretization_y = "<< CONFIG_sizeDiscretization_y << endl;

	cout << "OPTIMIZATION_nonuniform_t_sinh_scale_param  = "<< OPTIMIZATION_nonuniform_t_sinh_scale_param  << endl;
	cout << "OPTIMIZATION_nonuniform_x_sinh_scale_param  = "<< OPTIMIZATION_nonuniform_x_sinh_scale_param  << endl;
	cout << "OPTIMIZATION_nonuniform_y_sinh_scale_param  = "<< OPTIMIZATION_nonuniform_y_sinh_scale_param  << endl;


	//! Feller condition 
	double feller_condition = 2*Heston_params_kappa* Heston_params_theta - Heston_params_sigma*Heston_params_sigma; 
	if(feller_condition > 0)
	{
		cout << "Feller condition satisfied: feller = " << feller_condition << endl;
	}
	else
	{
		cout << "Feller condition NOT satisfied: feller = " << feller_condition << endl;
	}
	

	Interpolator_D2 interp_1 = test_Heston_PDE_solution(50);
	Interpolator_D2 interp_2 = test_Heston_PDE_solution(100);
	Interpolator_D2 interp_3 = test_Heston_PDE_solution(200);
	Interpolator_D2 interp_4 = test_Heston_PDE_solution(400);
      
	//interp_1.grid_value_.print();

	int N_x = 100;
	int N_y = 250;

	double x_start = 0.00001;
	double x_end   = 199.0;
	double x_step = (x_end-x_start)/(N_x-1);
	
	double y_start = 0.000000001;
	double y_end   = 0.02;
	double y_step = (y_end-y_start)/(N_y-1);

	vector<double> x_grid(N_x);
	x_grid[0] = x_start;
	for(int i=1; i<N_x; ++i)
	{
		x_grid[i] = x_grid[i-1] + x_step;
	}

	vector<double> y_grid(N_y);
	y_grid[0] = y_start;
	for(int i=1; i<N_y; ++i)
	{
		y_grid[i] = y_grid[i-1] + y_step;
	}
	
	Matrix result_1 = interp_1.interpolate(x_grid, y_grid);
    Matrix result_2 = interp_2.interpolate(x_grid, y_grid);
	Matrix result_3 = interp_3.interpolate(x_grid, y_grid);
	Matrix result_4 = interp_4.interpolate(x_grid, y_grid);
	Matrix result_1_2   = result_2 - result_1;
	Matrix result_2_3   = result_3 - result_2;
	Matrix result_3_4   = result_4 - result_3;

	std::string PDE_error_output_file_name_12 = DEBUG_output_path + "compare_PDE_result_12.csv";
	std::string PDE_error_output_file_name_23 = DEBUG_output_path + "compare_PDE_result_23.csv";
	std::string PDE_error_output_file_name_34 = DEBUG_output_path + "compare_PDE_result_34.csv";
	print_result(PDE_error_output_file_name_12 , x_grid, y_grid, result_1_2);
	print_result(PDE_error_output_file_name_23 , x_grid, y_grid, result_2_3);
	print_result(PDE_error_output_file_name_34 , x_grid, y_grid, result_3_4);
}

