//#include "Interpolator.h"
//#include "Params.h"
//#include <PDE/Matrix/Matrix.h>
//#include "useful_function.h"
//
//
//double f(double x, double y)
//{
//    return sin(2*x + 3.48*y);
//}
//
//void test_2D_Interpolator()
//{
//	//! grid
//	int N_x = 100;
//	int N_y = 50;
//
//	double x_start = -2;
//	double y_start = -3;
//
//	double x_step = 2*fabs(x_start)/(N_x-1);
//	double y_step = 2*fabs(y_start)/(N_y-1);
//
//	vector<double> x(N_x, 0.0); 
//	vector<double> y(N_y, 0.0); 
//	x[0] = x_start;
//    for(int i=1; i<N_x; ++i)
//	{
//	    x[i] = x[i-1] + x_step;
//	}
//	y[0] = y_start;
//    for(int i=1; i<N_y; ++i)
//	{
//	    y[i] = y[i-1] + y_step;
//	}
//
//	Matrix value(x.size(), y.size(), 0.0);
//	for(int i=0; i<x.size(); ++i)
//	{
//		for(int j=0; j<y.size(); ++j)
//		{
//		    value[i][j] = f(x[i],y[j]);
//		}
//	}
//
//	//! to interpolate 
//	int N_x_interp = N_x;
//	int N_y_interp = N_y;
//
//	double x_start_interp = x_start*0.87;
//	double y_start_interp = y_start*0.87;
//
//	double x_step_interp = 2*fabs(x_start_interp)/(N_x-1);
//	double y_step_interp = 2*fabs(y_start_interp)/(N_y-1);
//
//	vector<double> x_interp(N_x_interp, 0.0); 
//	vector<double> y_interp(N_y_interp, 0.0); 
//	x_interp[0] = x_start_interp;
//    for(int i=1; i<N_x_interp; ++i)
//	{
//	    x_interp[i] = x_interp[i-1] + x_step_interp;
//	}
//	y_interp[0] = y_start_interp;
//    for(int i=1; i<N_y_interp; ++i)
//	{
//	    y_interp[i] = y_interp[i-1] + y_step_interp;
//	}
//
//	//Matrix value_interp(x_interp.size(), y_interp.size(), 0.0);
//	
//
//	//! interpolator
//	Interpolator_D2 interp_2d(x,y,value);
//	Matrix value_interp = interp_2d.interpolate(x_interp,y_interp);
//
//	//! check result 
//	for(int i=0; i<x_interp.size(); ++i)
//	{
//		for(int j=0; j<y_interp.size(); ++j)
//		{
//		     value_interp[i][j] -= f(x_interp[i],y_interp[j]);
//		}
//
//	}
//	std::string PDE_error_output_file_name_1 = DEBUG_output_path + "interpolator_test_original.csv";
//	std::string PDE_error_output_file_name_2 = DEBUG_output_path + "interpolator_test.csv";
//	print_result(PDE_error_output_file_name_1 , x, y, value);
//	print_result(PDE_error_output_file_name_2 , x_interp, y_interp, value_interp);
//
//
//    
//}
//
