#include "Discretization.h"
#include <algorithm>
#include "useful_function.h"
#include "Range.h"
//#include <vector>


//const double epsilon_zero = 0.00000000001;
Discretization :: Discretization(  const Range& range_t,   int sizeDiscretization_t,
								   const Range& range_x,   int sizeDiscretization_x, 
								   const Range& range_y,   int sizeDiscretization_y, 
								   const std::vector<bool>& if_nonuniform_t, const std::vector<double>& t_nonuniform_center_c, const std::vector<double>& t_nonuniform_scale_param_c,
								   const std::vector<bool>& if_shift_t,		 const std::vector<double>& t_shifting_center, 
								   const std::vector<bool>& if_nonuniform_x, const std::vector<double>& x_nonuniform_center_c, const std::vector<double>& x_nonuniform_scale_param_c,
								   const std::vector<bool>& if_shift_x,		 const std::vector<double>& x_shifting_center, 
								   const std::vector<bool>& if_nonuniform_y, const std::vector<double>& y_nonuniform_center_c, const std::vector<double>& y_nonuniform_scale_param_c,
								   const std::vector<bool>& if_shift_y,		 const std::vector<double>& y_shifting_center
								 )  
								 :
								   discret_t(range_t,sizeDiscretization_t, 
								             if_nonuniform_t, t_nonuniform_center_c, t_nonuniform_scale_param_c,
								             if_shift_t, t_shifting_center),

								   discret_x(range_x,sizeDiscretization_x, 
								             if_nonuniform_x, x_nonuniform_center_c, x_nonuniform_scale_param_c,
								             if_shift_x, x_shifting_center),

								   discret_y(range_y,sizeDiscretization_y,
											 if_nonuniform_y, y_nonuniform_center_c, y_nonuniform_scale_param_c,
								             if_shift_y, y_shifting_center)
{}


//discretization discretization :: create_new_discretization_by_refine(int refine_ratio_time, int refine_ratio_s, int refine_ratio_v)   //! refine each interval to refine_ration times, equal-distance;
//{
//    throw("Error in function discretization discretization :: create_new_discretization_by_refine(...), function not implemented.");
//    /*return discretization(  range_t_,   (sizeDiscretization_t_-1)*refine_ratio_time  + 1,
//						    range_x_,   (sizeDiscretization_x_-1)*refine_ratio_space + 1, 
//							if_uniform_x_, x_nonuniform_center_c_ , x_nonuniform_scale_param_c_,
//							if_shift_x_,   x_shifting_center_);*/
//}

Discretization Discretization::refine(int ratio_t, int ratio_x, int ratio_y) const
{
	if(ratio_t<=0 || ratio_x<=0 || ratio_y<=0)
	{
	    throw ("Error in function discretization::refine(...), invalid ratio value.");
	}

	int sizeDiscretization_t_2 = (discret_t.sizeDiscretization_-1)*ratio_t + 1; 
	int sizeDiscretization_x_2 = (discret_x.sizeDiscretization_-1)*ratio_x + 1; 
	int sizeDiscretization_y_2 = (discret_y.sizeDiscretization_-1)*ratio_y + 1;

	return Discretization(  discret_t.get_range(),  sizeDiscretization_t_2,
						    discret_x.get_range(),  sizeDiscretization_x_2, 
						    discret_y.get_range(),  sizeDiscretization_y_2,

							discret_t.if_nonuniform_, discret_t.nonuniform_center_c_, discret_t.nonuniform_scale_param_c_,
							discret_t.if_shift_,      discret_t.shifting_center_, 

						    discret_x.if_nonuniform_, discret_x.nonuniform_center_c_, discret_x.nonuniform_scale_param_c_,
						    discret_x.if_shift_,      discret_x.shifting_center_, 

						    discret_y.if_nonuniform_, discret_y.nonuniform_center_c_, discret_y.nonuniform_scale_param_c_,
						    discret_y.if_shift_,      discret_y.shifting_center_
						  );
}


void Discretization :: print() const
{
    discret_t.print();
	discret_x.print();
	discret_y.print();
}

void Discretization ::print_t() const
{
	discret_t.print();
}

void Discretization ::print_x() const
{
	discret_x.print();
}

void Discretization ::print_y() const
{
	discret_y.print();
}