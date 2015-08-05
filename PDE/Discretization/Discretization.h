#pragma once

#include <iostream>
//#include <PDE/Discretization/Range.h>
#include "Discretization_Simple.h"
#include <vector>
#include <boost/shared_ptr.hpp>

//! TODO: too many reboundant code, need to restructure the Disretization class as one class contains 3 attribute, 
//!       each contains 1D discretization: time, space, volatility.

//! We suppose the schemas is time homogenerous.
//! If not, it will cause problem in "Schemas.h", in function "initializeSchemas(double t_n, double t_nPlus)",
//! when call: "initlialize_first_last_column_elements(t_n, t_nPlus, t_x_first, x_last" where "t_x_first/last"
//! is time independent ! 


//! X_tilde is a part of X (without first and last rows)

//! Suppose always that: the domain of discretization is rectangular
class Discretization
{
public:
	Discretization_Simple discret_t;
	Discretization_Simple discret_x;
	Discretization_Simple discret_y;

public:
	//! CONSTRUCTOR
	//! uniform discretization

	Discretization(const Range& range_t,   int sizeDiscretization_t,
				   const Range& range_x,   int sizeDiscretization_x, 
				   const Range& range_y,   int sizeDiscretization_y,

				   const std::vector<bool>& if_nonuniform_t = std::vector<bool>(1,false), const std::vector<double>& t_nonuniform_center_c = std::vector<double>(1,-99999.9), const std::vector<double>& t_nonuniform_scale_param_c = std::vector<double>(1,-99999.9),
				   const std::vector<bool>& if_shift_t      = std::vector<bool>(1,false), const std::vector<double>& t_shifting_center     = std::vector<double>(1,-99999.9), 

				   //! non-uniform / shifting paramters for s
				   const std::vector<bool>& if_nonuniform_x = std::vector<bool>(1,false), const std::vector<double>& x_nonuniform_center_c = std::vector<double>(1,-99999.9), const std::vector<double>& x_nonuniform_scale_param_c = std::vector<double>(1,-99999.9),
				   const std::vector<bool>& if_shift_x      = std::vector<bool>(1,false), const std::vector<double>& x_shifting_center     = std::vector<double>(1,-99999.9), // shifting_center = K;

				   //! non-uniform / shifting paramters for v				   
				   const std::vector<bool>& if_nonuniform_y = std::vector<bool>(1,false), const std::vector<double>& y_nonuniform_center_c = std::vector<double>(1,-99999.9), const std::vector<double>& y_nonuniform_scale_param_c = std::vector<double>(1,-99999.9),
				   const std::vector<bool>& if_shift_y      = std::vector<bool>(1,false), const std::vector<double>& y_shifting_center     = std::vector<double>(1,-99999.9) 
				   ); 

	//discretization create_new_discretization_by_refine(int refine_ratio_time, int refine_ratio_x, int refine_ratio_y); //! refine each interval to refine_ration times, equal-distance;

	//! DESTRUCTOR
	virtual ~Discretization(){};  //! no pointer, so need of destructor ...


	std::vector<double> get_discret_t() const {return discret_t.get_discret();}
	std::vector<double> get_discret_x() const {return discret_x.get_discret();}
	std::vector<double> get_discret_y() const {return discret_y.get_discret();}

	std::vector<double> get_discret_x_tilde() const {return discret_x.get_discret_tilde();}
	std::vector<double> get_discret_y_tilde() const {return discret_y.get_discret_tilde();}

	double get_discret_t       (int index) const {return discret_t.get_discret(index);}
	double get_discret_x       (int index) const {return discret_x.get_discret(index);}
	double get_discret_x_tilde (int index) const {return discret_x.get_discret_tilde(index);}
	double get_discret_y       (int index) const {return discret_y.get_discret(index);}
	double get_discret_y_tilde (int index) const {return discret_y.get_discret_tilde(index);}

	//! parameter t, not used here!
	double get_discret_x_left_border (int t_index)  const {return discret_x.get_discret_left_border(t_index);}
	double get_discret_x_right_border(int t_index)  const {return discret_x.get_discret_right_border(t_index);}
	double get_discret_y_left_border (int t_index)  const {return discret_y.get_discret_left_border(t_index);}
	double get_discret_y_right_border(int t_index)  const {return discret_y.get_discret_right_border(t_index);}

	int get_sizeDiscret_t()        const {return (int)discret_t.get_sizeDiscret();}
	int get_sizeDiscret_x()        const {return (int)discret_x.get_sizeDiscret();}
	int get_sizeDiscret_x_tilde()  const {return (int)discret_x.get_sizeDiscret_tilde();}
	int get_sizeDiscret_y()        const {return (int)discret_y.get_sizeDiscret();}
	int get_sizeDiscret_y_tilde()  const {return (int)discret_y.get_sizeDiscret_tilde();}

	double get_delta_t(int index) const {return discret_t.get_discret_delta(index);}
	double get_delta_x(int index) const {return discret_x.get_discret_delta(index);}
	double get_delta_y(int index) const {return discret_y.get_discret_delta(index);}

	//! only for uniform discretization ...
	double get_discret_t_pas() const{return discret_t.get_discret_pas();}
	double get_discret_x_pas() const{return discret_x.get_discret_pas();}
	double get_discret_y_pas() const{return discret_y.get_discret_pas();}

	//! get simple discretization reference
	Discretization_Simple& get_t_Discretization_Simple_ref(){return discret_t;}
	Discretization_Simple& get_x_Discretization_Simple_ref(){return discret_x;}
	Discretization_Simple& get_y_Discretization_Simple_ref(){return discret_y;}

	//! for extrapolation: grid refinement 
    Discretization refine(int ratio_t, int ratio_x, int ratio_y) const;

	void print() const;
	void print_t() const;
	void print_x() const;
	void print_y() const;
};

typedef boost::shared_ptr<Discretization> Discretization_PTR;
typedef boost::shared_ptr<const Discretization> Discretization_CONSTPTR;