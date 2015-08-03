////! Fwd PDE
//
//#include <iostream>
//#include <fstream>
//#include <fstream>
//#include <string>
//#include <sstream>
//#include <ctime>
//#include <cmath>
//
//#include "Discretization.h"
//#include "BoundaryCondition_Heston_Put.h"
//#include "BoundaryCondition_D2_density.h"
//
//#include "PDE_ADI_Solver.h"
////#include "PDE_ADI_Extrapolation_Solver.h"
//
//#include "Scheme_ADI.h"
//#include "Scheme_Yanenko.h"  // 1
//#include "Scheme_Douglas.h"  // 2 
//#include "Scheme_CS.h"       // 3
//#include "Scheme_MCS.h"      // 4
//#include "Scheme_HV.h"       // 5
//
////#include "HestonModel.h"
////#include "Model_discret.h"
////#include "Model.h"
//#include "Model_discret.h"
//#include "Heston_analytical_pricer.h"
//#include "useful_function.h"
//
//#include <PDE/Matrix/Matrix.h>
//#include "Params.h"
//#include <fstream>
//
//#include "Interpolator.h"
///*
//Heston Model: 
//  \frac{dS}{S} = r*dt + sqrt(v)*dW^1
//  dv = \kappa(\theta - v)*dt + sigma*sqrt(v)*dW^2
//  <dW^1,dW^2> = \rho*dt
//*/
//
////! some constat about the stability: 
////! 1. When the sinh non-unifor scheme is used, n_t_discret = 100; n_x_discre = 200, n_y_discret = 50:
////! Yanenko is stable, Douglas (1) is stable,  CS (1, 1/2) is stable
////! Douglas(1/2) not stable !  CS (1/2, 1/2) not stable ! 
//
//
//const int    scheme_type           = 3;   // 1: Yanenko, 2:Douglas, 3:CS, 4:MCS, 5:HV
//const double scheme_douglas_theta  = 0.5; // used for scheme: 2,3,4,5
//const double scheme_douglas_lambda = 0.5; // used for scheme: 3,4
//
//
//
//const bool  DEBUG_test_PDE_solver_extrapolation   = false; // Don't work ... ! What to do for 2D extrapolation ???
//
//const bool  DEBUG_PDE_test_only_PDE				  = true; 
//const bool  DEBUG_PDE_test_PDE_vs_analytical	  = !DEBUG_PDE_test_only_PDE;
//
//const bool  DEBUG_PDE_error_print_to_file		  = true;  //! for all grid (S,v)
//const bool  DEBUG_PDE_error_print_to_screen		  = false;  //! for only one grid :)
//
//const double CONFIG_S_max = fwd_PDE_SLV_SABR_S0*10;//Option_params_K*8.0;
//const double CONFIG_v_max = fwd_PDE_SLV_SABR_v0*5;//5.0;
//
//const int    CONFIG_multiplicator = 1;
////! you should normally use config: t=50, x=100, y =50
//const int    CONFIG_sizeDiscretization_t = 50*CONFIG_multiplicator;
//const int    CONFIG_sizeDiscretization_x = 100*CONFIG_multiplicator; //Attention: when compare the precision of uniform and non-uniform result, we need a grid at S = 100 !
//const int    CONFIG_sizeDiscretization_y = 50*CONFIG_multiplicator; 
//
//const bool   OPTIMIZATION_if_nonuniform_sinh_grid_t = false;   
//const bool   OPTIMIZATION_if_shfit_grid_t           = false;   
//const bool   OPTIMIZATION_if_nonuniform_sinh_grid_x = false;   // nonuniform discretization on space: x-direction
//const bool   OPTIMIZATION_if_shfit_grid_x           = false;   // OPTIMIZATION_if_nonuniform_sinh_grid_x;
//const bool   OPTIMIZATION_if_nonuniform_sinh_grid_y = false;   // 0onuniform discretization on space: y-direction
//const bool   OPTIMIZATION_if_shfit_grid_y           = false;
//
//const double OPTIMIZATION_nonuniform_t_sinh_scale_param = 1/10000.0;
////const double OPTIMIZATION_nonuniform_x_sinh_scale_param = Option_params_K/5; // Heston pricing
//const double OPTIMIZATION_nonuniform_x_sinh_scale_param = 1/100.0; // Fwd PDE
//const double OPTIMIZATION_nonuniform_y_sinh_scale_param = 1/100.0;
//
//
//Interpolator_D2 test_SLV_PDE_solution(int CONFIG_sizeDiscretization_y)
//{	
//	std::cout << "Model and Product Direction_Parameters" << std::endl;
//
//	double r = SLV_SABR_params_r; 
//
//	//! Heston Direction_Parameters
//	double alpha = SLV_SABR_params_alpha;
//	double beta  = SLV_SABR_params_beta;
//	double rho   = SLV_SABR_params_rho;
//
//	//! Product Direction_Parameters
//	double K = Option_params_K;
//	double T = Option_params_T; 
//
//	//! FD param
//	double S_min = 0.0;           // S_max = 2*K is not efficient! 
//	double S_max = CONFIG_S_max;  // S_max = 2*K is not efficient! 
//	double v_min = 0.0;
//	double v_max = CONFIG_v_max;
//
//
//	//! resume: the biggest error is the time! 
//	int sizeDiscretization_t = CONFIG_sizeDiscretization_t + 1;  
//	int sizeDiscretization_x = CONFIG_sizeDiscretization_x + 1; 
//	int sizeDiscretization_y = CONFIG_sizeDiscretization_y + 1;
//
//	//! Range and discretization
//	Range range_t(0.0,  T);        
//	Range range_x(S_min,  S_max);   
//	Range range_y(v_min,  v_max);  
//
//	//! non-uniform 
//		//! direction t
//	int t_discretization_size = 1;
//	std::vector<bool>   if_nonuniform_t(t_discretization_size);//              = OPTIMIZATION_if_nonuniform_sinh_grid_x;
//    std::vector<bool>   if_shift_t(t_discretization_size); //                  = OPTIMIZATION_if_shfit_grid_x;
//	std::vector<double> t_nonuniform_center_c(t_discretization_size);//        = K; 
//	std::vector<double> t_nonuniform_scale_param_c(t_discretization_size);//   = K/5.0; 
//	std::vector<double> t_shifting_center(t_discretization_size); //            = K;
//
//	if_nonuniform_t[0]              = OPTIMIZATION_if_nonuniform_sinh_grid_t;
//    if_shift_t[0]                   = OPTIMIZATION_if_shfit_grid_t;
//	t_nonuniform_center_c[0]        = 0.0; 
//	t_nonuniform_scale_param_c[0]   = OPTIMIZATION_nonuniform_t_sinh_scale_param; //1.0/500.0; 
//	t_shifting_center[0]            = -999999.9;
//
//	if(if_fwd_PDE == true)
//	{
//		////! add second time
//  //     	if_nonuniform_t.push_back(false);             //           = true;//OPTIMIZATION_if_nonuniform_sinh_grid_t;
//		//if_shift_t.push_back(false)    ;              // = OPTIMIZATION_if_shfit_grid_t;
//		//t_nonuniform_center_c.push_back(0.0);		  //        = 0.0; 
//		//t_nonuniform_scale_param_c.push_back(1.0/50); //   = 1.0/25.0; 
//		//t_shifting_center.push_back(-99999.9);        //            = -999999.9;
//	}
//
//		//! direction x
//	int x_discretization_size = 1;
//	std::vector<bool>   if_nonuniform_x(x_discretization_size);//              = OPTIMIZATION_if_nonuniform_sinh_grid_x;
//    std::vector<bool>   if_shift_x(x_discretization_size); //                  = OPTIMIZATION_if_shfit_grid_x;
//	std::vector<double> x_nonuniform_center_c(x_discretization_size);//        = K; 
//	std::vector<double> x_nonuniform_scale_param_c(x_discretization_size);//   = K/5.0; 
//	std::vector<double> x_shifting_center(x_discretization_size); //            = K;
//
//	if_nonuniform_x[0]              = OPTIMIZATION_if_nonuniform_sinh_grid_x;
//    if_shift_x[0]                   = OPTIMIZATION_if_shfit_grid_x;
//	x_nonuniform_center_c[0]        = K; 
//	x_nonuniform_scale_param_c[0]   = K/5.0; 
//	x_shifting_center[0]            = K;
//
//	if(if_fwd_PDE == true)
//	{
//		x_nonuniform_center_c[0]      = fwd_PDE_S0; 
//		x_nonuniform_scale_param_c[0] = OPTIMIZATION_nonuniform_x_sinh_scale_param; 
//		x_shifting_center[0]          = fwd_PDE_S0;
//	}
//
//		//! direction y
//	int y_discretization_size = 1;
//	std::vector<bool>   if_nonuniform_y(y_discretization_size);		    
//    std::vector<bool>   if_shift_y(y_discretization_size);				
//	std::vector<double> y_nonuniform_center_c(y_discretization_size);   
//	std::vector<double> y_nonuniform_scale_param_c(y_discretization_size); 
//	std::vector<double> y_shifting_center(y_discretization_size);          
//
//	if_nonuniform_y[0]		      = OPTIMIZATION_if_nonuniform_sinh_grid_y;
//    if_shift_y[0]				  = OPTIMIZATION_if_shfit_grid_y; //! for Fwd PDE
//	y_nonuniform_center_c[0]      = 0.0;     
//	y_nonuniform_scale_param_c[0] = OPTIMIZATION_nonuniform_y_sinh_scale_param;  // V_max/500 is the experienced choice! 
//	y_shifting_center[0]          = -99999.9;
//
//	if(if_fwd_PDE == true)
//	{
//		//! first dimension: 
//		if_nonuniform_y[0]			  = OPTIMIZATION_if_nonuniform_sinh_grid_y; // OPTIMIZATION_if_nonuniform_sinh_grid_y;
//		if_shift_y[0]				  = OPTIMIZATION_if_shfit_grid_y; // OPTIMIZATION_if_shfit_grid_y; //! for Fwd PDE
//		y_nonuniform_center_c[0]      = 0.0;//Heston_params_theta;     
//		y_nonuniform_scale_param_c[0] = OPTIMIZATION_nonuniform_y_sinh_scale_param; // OPTIMIZATION_nonuniform_y_sinh_scale_param; //1.0/10000.0;  // V_max/500 is the experienced choice! 
//		y_shifting_center[0]          = fwd_PDE_v0;
//
//		////! second dimension:
//		//if_nonuniform_y.push_back(false);
//		//if_shift_y.push_back(true);					                  
//		//y_nonuniform_center_c.push_back(fwd_PDE_v0);			      
//		//y_nonuniform_scale_param_c.push_back(1.0/*OPTIMIZATION_nonuniform_y_sinh_scale_param*/);
//		//y_shifting_center.push_back(fwd_PDE_v0); 
//
//		////! third
//		//if_nonuniform_y.push_back(false);
//		//if_shift_y.push_back(true);					                  
//		//y_nonuniform_center_c.push_back(Heston_params_theta);			      
//		//y_nonuniform_scale_param_c.push_back(1.0/5000.0);
//		//y_shifting_center.push_back(fwd_PDE_v0);
//	}
//
//	Discretization discret( range_t, sizeDiscretization_t,  
//							range_x, sizeDiscretization_x,  
//							range_y, sizeDiscretization_y,
//							if_nonuniform_t,  t_nonuniform_center_c,  t_nonuniform_scale_param_c,
//							if_shift_t,    t_shifting_center,
//							if_nonuniform_x,  x_nonuniform_center_c,  x_nonuniform_scale_param_c,
//							if_shift_x,    x_shifting_center,
//							if_nonuniform_y,  y_nonuniform_center_c,  y_nonuniform_scale_param_c,
//							if_shift_y,    y_shifting_center);  
//
//	//! SABR Model
//	SABRModel model  (alpha,beta,rho,r);
//	//SABR_Model_discret sabr_model_discret(model);
//	//! SABR discret Model 
//	SLV_SABR_Model_discret slv_sabr_model_discret(model,discret);
//
//
//    //! Term Structure, Payoff 	&&  Boundary Condition
//	Term_Structure* cts  = new Const_Term_Structure(r);
//	PayoffYuan* payoff   = new Put(K);
//
//	//! Backward PDE
//	//BoundaryCondition_Heston_Put bc(cts, payoff, model, discret);
//	//! Forward PDE
//	BoundaryCondition_D2_density bc(cts, payoff, model, discret ,fwd_PDE_SLV_SABR_S0, fwd_PDE_SLV_SABR_v0);
//
//	//! --------------------------------------------------
//	//! 
//    //!           Test PDE solver 
//	//! 
//	//! ---------------------------------------------------
//	//! new Scheme class 
//	Scheme_ADI* scheme_new_ptr = NULL;
//	if(scheme_type ==1)
//	{
//		scheme_new_ptr= new Scheme_Yanenko(
//			discret,
//			//sabr_model_discret,
//			slv_sabr_model_discret,
//			bc
//			);
//	}
//	else if(scheme_type ==2)
//	{
//		scheme_new_ptr = new Scheme_Douglas(
//			discret,
//			//sabr_model_discret,
//			slv_sabr_model_discret,
//			bc,
//			scheme_douglas_theta
//			);
//	}
//	else if (scheme_type == 3)
//	{
//		scheme_new_ptr = new Scheme_CS(
//			discret,
//			//sabr_model_discret,
//			slv_sabr_model_discret,
//			bc,
//			scheme_douglas_theta,
//			scheme_douglas_lambda
//			);
//	}
//	else if(scheme_type == 4)
//	{
//		scheme_new_ptr = new Scheme_MCS(
//			discret,
//			//sabr_model_discret,
//			slv_sabr_model_discret,
//			bc,
//			scheme_douglas_theta,
//			scheme_douglas_lambda
//			);
//	}
//	else if(scheme_type == 5)
//	{
//		scheme_new_ptr = new Scheme_HV(
//			discret,
//			//sabr_model_discret,
//			slv_sabr_model_discret,
//			bc,
//			scheme_douglas_theta
//			);
//	}
//	else
//	{
//	    throw ("Error in Test_PDE.cpp, scheme_type  does not take valid value");
//	}
//
//	PDE_ADI_Solver   pde_solver(*scheme_new_ptr);
//
//	clock_t t1 = clock();
//	pde_solver.solve_PDE();	
//	clock_t t2 = clock();
//	cout << "time = " << (double)(t2-t1)/CLOCKS_PER_SEC << " seconds"<< endl;
//
//	Matrix Fwd_PDE_result (pde_solver.get_result());
//    
//	Matrix Fwd_PDE_result_x(discret.get_sizeDiscret_x_tilde(),3,0.0); // totoal proba, positive proba, negative proba ! 
//
//	double sum_proba          = 0.0;
//	double sum_proba_positive = 0.0;
//	double sum_proba_negative = 0.0;
//	for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
//	{
//		double sum_proba_x = 0.0;
//		double sum_proba_positive_x = 0.0;
//		double sum_proba_negative_x = 0.0;
//
//		int index_i = i+1;
//		for( int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
//		{
//			int index_j = j+1;
//            double density = Fwd_PDE_result[index_i][index_j]; // then it is the density
//			double coeff   = 1.0/4
//		                    *( discret.get_delta_x(i)+discret.get_delta_x(i+1) )
//						    *( discret.get_delta_y(j)+discret.get_delta_y(j+1) );
//
//			double coeff_y = 1.0/2*( discret.get_delta_y(j)+discret.get_delta_y(j+1) );
//
//			double proba   = coeff*density;  
//			if(density>0.0)
//			{
//				sum_proba_positive += coeff*density;  
//				//! x
//				sum_proba_positive_x += coeff_y*density;
//			}
//			else
//			{
//			    sum_proba_negative  += coeff*density;  
//				//! x
//				sum_proba_negative_x  += coeff_y*density;  
//			}
//		}
//		sum_proba_x = sum_proba_positive_x + sum_proba_negative_x;
//
//		//! save backup 
//		Fwd_PDE_result_x[i][0] = sum_proba_x;
//		Fwd_PDE_result_x[i][1] = sum_proba_positive_x;
//		Fwd_PDE_result_x[i][2] = sum_proba_negative_x;
//	}
//	sum_proba = sum_proba_positive + sum_proba_negative;
//
//	cout << " ---- |||| ---- |||| Sum_Proba          = " << sum_proba << endl;  //<< ",  sum_proba_positive = " << sum_proba_positive << " + " << "sum_proba_negative = " << sum_proba_negative << endl;
//	cout << " ---- |||| ---- |||| Sum_Proba_positive = " << sum_proba_positive << endl;
//	cout << " ---- |||| ---- |||| Sum_Proba_negative = " << sum_proba_negative << endl;
//	cout << " ---- |||| ---- |||| Sum_total_error    = " << fabs(sum_proba_positive - 1) + fabs(sum_proba_negative) << endl;
//
//
//
//	////! analytical price
//	//int num_K = 17;
//	//double step = 0.05;
//	////double K_start = Option_params_K*(num_K-1)/
//	//vector<double> v_K(num_K);
//	//vector<double> v_analytical_price(num_K);
//
//	//cout << "strike vector: "<< endl;
//	//cout << Option_params_K << endl;
//	//v_K[0] = Option_params_K;
//	//for(unsigned int i=1; i<v_K.size()-1; i+=2)
//	//{
//	//	v_K[i]   =  Option_params_K *(1+i*step);
//	//	v_K[i+1] =  Option_params_K *(1-i*step);
//	//	cout << v_K[i] << endl;
//	//	cout << v_K[i+1] << endl;
//	//}
//
//	//for(unsigned int i=0; i<v_K.size(); ++i)
//	//{
//	//	Heston_analytical_pricer pricer_ana(v_K[i]);
//	//	v_analytical_price[i] = pricer_ana.heston_put(fwd_PDE_S0,fwd_PDE_v0, 1024*10);
//
//	//}
//
//	//vector<double> v_numerical_price(num_K);
//	////! first numerical integration price
//	//for(unsigned int kk=0; kk<v_K.size()-1; ++kk)
//	//{
//	//	    double KK = v_K[kk];
//	//		double sum_price          = 0.0;
//	//		double sum_price_positive = 0.0;
//	//		double sum_price_negative = 0.0;
//	//		for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
//	//		{
//	//			int index_i = i+1;
//	//			double payoff = 1;
//	//			double x = discret.get_discret_x(index_i);
//	//			if(x<KK)
//	//				payoff = KK - x;
//	//			else
//	//				payoff = 0.0;
//
//	//			for( int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
//	//			{
//	//				int index_j = j+1;
//	//				double density = Fwd_PDE_result[index_i][index_j]; // then it is the density
//	//				double coeff   = 1.0/4
//	//								*( discret.get_delta_x(i)+discret.get_delta_x(i+1) )
//	//								*( discret.get_delta_y(j)+discret.get_delta_y(j+1) );
//
//	//				double proba   = coeff*density;  
//	//				if(density>0.0)
//	//				{
//	//					sum_price_positive += coeff*density*payoff;  
//	//				}
//	//				else
//	//				{
//	//				   sum_price_negative  += coeff*density*payoff;  
//	//				}
//	//			}
//	//		}
//	//		sum_price_positive *= exp(-Option_params_T*Heston_params_r);
//	//		sum_price_negative *= exp(-Option_params_T*Heston_params_r);
//	//		sum_price = sum_price_positive + sum_price_negative;
//	//		v_numerical_price[kk] = sum_price;
//	//		cout << "1st numerical price = " << sum_price << endl;
//	//		//cout << "price_numerical_positive  = " << sum_price_positive << endl;
//	//		//cout << "price_numerical_negative  = " << sum_price_negative << endl;
//
//
//	//		//!second numerical integration price
//	//		vector<double> x_grid(discret.get_discret_x_tilde());
//	//		vector<double> y_grid(discret.get_discret_y_tilde());
//
//	//		vector<double> x_value(discret.get_sizeDiscret_x_tilde(),0.0);
//	//		vector<double> y_value(discret.get_sizeDiscret_y_tilde(),0.0);
//
//	//		for(int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
//	//		{
//	//			int index_i = i+1;
//	//			double payoff = 1;
//	//			//double x = discret.get_discret_x(index_i);
//	//			//if(x<KK)
//	//			//	payoff = KK - x;
//	//			//else
//	//			//	payoff = 0.0;
//
//	//			for(int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
//	//			{
//	//				 int index_j = j+1;
//	//				 y_value[j] = Fwd_PDE_result[index_i][index_j];
//	//			}
//	//			Interpolator interp_y(y_grid,y_value);
//	//			x_value[i] = interp_y.calculate_integral_approximation()*payoff;
//	//		}
//	//		Interpolator interp_x(x_grid,x_value);
//	//		double second_analytical_price = interp_x.calculate_integral_approximation();
//	//		//second_analytical_price *= exp(-Option_params_T*Heston_params_r);
//	//		cout << "2nd numerical price = " << second_analytical_price << endl;
//	//		
//
//	//		//! third numerical integration price
//	//		double third_sum_price          = 0.0;
//	//		vector<double> v_sum_x(discret.get_sizeDiscret_x_tilde(),0.0);
//	//		vector<double> v_sum_x_positive(discret.get_sizeDiscret_x_tilde(),0.0);
//	//		vector<double> v_sum_x_negative(discret.get_sizeDiscret_x_tilde(),0.0);
//	//		for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
//	//		{
//	//			int index_i = i+1;
//	//			double payoff = 1;
//	//			//double x = discret.get_discret_x(index_i);
//	//			//if(x<KK)
//	//			//	payoff = KK - x;
//	//			//else
//	//			//	payoff = 0.0;
//	//			double sum_x          = 0.0;
//	//			double sum_x_positive = 0.0;
//	//			double sum_x_negative = 0.0;
//	//			for( int j=0; j<discret.get_sizeDiscret_y_tilde(); ++j)
//	//			{
//	//				int index_j = j+1;
//	//				double density = Fwd_PDE_result[index_i][index_j]; // then it is the density
//	//				double coeff   = ( discret.get_delta_y(j)+discret.get_delta_y(j+1) )/2;
//
//	//				double proba   = coeff*density;  
//	//				if(density>0.0)
//	//				{
//	//					sum_x_positive += coeff*density*payoff;  
//	//				}
//	//				else
//	//				{
//	//				   sum_x_negative  += coeff*density*payoff;  
//	//				}
//	//				sum_x = sum_x_positive + sum_x_negative;
//	//			}
//	//			v_sum_x_positive[i] = sum_x_positive;
//	//			v_sum_x_negative[i] = sum_x_negative;
//	//			v_sum_x[i]          = sum_x;
//	//		}
//
//	//		double hehe = 0.0;
//	//		for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
//	//		{
//	//			double interval = ( discret.get_delta_x(i)+discret.get_delta_x(i+1) )/2.0;
//	//			hehe += v_sum_x[i]*interval;
//	//		}
//	//		cout << "3rd numerical price = " << hehe << endl;
//	//}
//	Interpolator_D2 interp_2d(discret.get_discret_x(), discret.get_discret_y(), Fwd_PDE_result);
//
//			////sum_price_positive *= exp(-Option_params_T*Heston_params_r);
//			////sum_price_negative *= exp(-Option_params_T*Heston_params_r);
//			//sum_price = sum_price_positive + sum_price_negative;
//			//cout << "price_numerical  = " << sum_price << endl;
//			//cout << "price_numerical_positive  = " << sum_price_positive << endl;
//			//cout << "price_numerical_negative  = " << sum_price_negative << endl;
//	
//
//	if (DEBUG_PDE_error_print_to_file == true)
//	{
//		std::stringstream ss1;
//		ss1 << CONFIG_sizeDiscretization_x;
//		string s1 = ss1.str();
//
//		std::stringstream ss2;
//		ss2 << CONFIG_sizeDiscretization_y;
//		string s2 = ss2.str();
//		std::string PDE_error_output_file_name = DEBUG_output_path + s1 + s2 + "PDE_error.csv";
//
//		Matrix PDE_result = Fwd_PDE_result;
//		//! get analytical result 
//		Heston_analytical_pricer analytical_pricer;
//
//		vector<double> discret_x(discret.get_discret_x());
//		vector<double> discret_y(discret.get_discret_y());
//
//
//		//! Fwd PDE
//		double Error_S_min = S_min;
//		double Error_S_max = S_max; // 1.5 * Option_params_K;
//		double Error_v_min = v_min;
//		double Error_v_max = v_max; //1.0;
//
//
//		vector<double> x_vector_value;  //! For analytical price
//		vector<int>    x_vector_index;  
//		vector<double> y_vector_value;  //! For PDE price (index) 
//		vector<int>    y_vector_index;   
//
//		for( int i=0; i<discret.get_sizeDiscret_x_tilde(); ++i)
//		{
//			double x = discret_x[i];
//		    if(x>=Error_S_min-0.000000001 && x <= Error_S_max+0.000000001)
//			{
//				x_vector_value.push_back(x);
//				x_vector_index.push_back(i);
//			}
//		}
//
//		for(unsigned int i=0; i<discret_y.size(); ++i)
//		{
//			double y = discret_y[i];
//		    if(y>=Error_v_min-0.000000001 && y <= Error_v_max+0.000000001)
//			{
//				y_vector_value.push_back(y);
//				y_vector_index.push_back(i);
//			}
//		}
//		
//		//Matrix analytical_price(1,1,0.0);
//		//if(DEBUG_PDE_test_PDE_vs_analytical == true)
//		//{
//		//	analytical_price  = analytical_pricer.heston_put_v(x_vector_value,y_vector_value);
//		//}
//		Matrix PDE_error((int)x_vector_value.size(), (int)y_vector_value.size(),-9999.0);  
//		for(unsigned int i=0; i<x_vector_value.size(); ++i)
//		{
//			for(unsigned int j=0; j<y_vector_value.size(); ++j)
//			{
//			    if(DEBUG_PDE_test_PDE_vs_analytical == true)
//				{
//				   // PDE_error[i][j] = PDE_result[x_vector_index[i]][y_vector_index[j]] - analytical_price [i][j];
//				}
//				else if (DEBUG_PDE_test_only_PDE == true)
//				{
//				    PDE_error[i][j] = PDE_result[i][j];
//				}
//			}
//		}
//		print_result(PDE_error_output_file_name , x_vector_value, y_vector_value, PDE_error);
//
//		std::string PDE_error_output_file_name_x = DEBUG_output_path + "to_compare_with_1D_localVol_PDE_result.csv";
//		print_result(PDE_error_output_file_name_x , x_vector_value, y_vector_value, Fwd_PDE_result_x);
//	}
//
//	//for(unsigned int i=0; i<v_K.size(); ++i)
//	//{
//	//	cout << "v_K[" << i << "] = "<< v_K[i] << " 's price: " << v_analytical_price[i] << " ~ " << v_numerical_price[i] << endl;
//	//}
//
//	delete cts;
//	delete payoff;
//	//delete scheme_ptr;
//	delete scheme_new_ptr;
//
//	return interp_2d;
//}
//
//
//void test_SLV_PDE()
//{	
//	Interpolator_D2 interp_1 = test_SLV_PDE_solution(CONFIG_sizeDiscretization_y);
//	//cout << "cumulated probability = " << interp_1.<< endl;
//}
//
