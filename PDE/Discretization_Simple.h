#pragma once

#include <iostream>
#include "Range.h"
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
class Discretization_Simple
{
public:
	Range  range_;
	int    sizeDiscretization_;
								  
	//! non-uniform
	std::vector<bool>   if_nonuniform_;
	std::vector<double> nonuniform_center_c_;
	std::vector<double> nonuniform_scale_param_c_;

	//! shifting 
	std::vector<bool>   if_shift_;
	std::vector<double> shifting_center_; 

	std::vector<double> discretization_;
	std::vector<double> discretization_tilde_;  // size = x.size()-2 (without the firste and the last)
	std::vector<double> discretization_delta_;  // discretization_x_delta_[i] = discretization_x_[i+1] - discretization_x_[i]
												// size = discretization_x_.size()-

public:

	//! for s non-uniform
	Discretization_Simple( const Range& range,    int sizeDiscretization,
						   //! non-uniform 
						   const std::vector<bool>& if_nonuniform = std::vector<bool>(1,false), const std::vector<double>& nonuniform_center_c = std::vector<double>(1,-99999.9), const std::vector<double>& nonuniform_scale_param_c = std::vector<double>(1,-99999.9),
						   //! shifting paramters
						   const std::vector<bool>& if_shift = std::vector<bool>(1,false),      const std::vector<double>& shifting_center = std::vector<double>(1,-99999.9)
						   ); 


	std::vector<double> discret(int sizeDiscretization, int index);
	//Discretization_Simple create_new_discretization_by_refine(int refine_ratio); //! refine each interval to refine_ration times, equal-distance;

	virtual ~Discretization_Simple(){};  //! no pointer, so need of destructor ...
	void initialize_other_tables();

	//! discretization uniform 
	void uniform_discretisize_vector(Range& range, std::vector<double>& dicret) const;


	void add_discret(std::vector<double> & discret_v); // not efficient, but I don't care :)
	std::vector<double> get_discret()       const {return discretization_;}
	std::vector<double> get_discret_tilde() const {return discretization_tilde_;}
	std::vector<double> get_discret_delta() const {return discretization_delta_;}

	double get_discret       (int index) const;
	double get_discret_tilde (int index) const;
	double get_discret_delta (int index) const; //{return discretization_delta_[index];}

	//! parameter t, not used here!
	double get_discret_left_border (int t_index)  const {return discretization_[0];}
	double get_discret_right_border(int t_index)  const {return discretization_[discretization_.size()-1];}

	int get_sizeDiscret()        const {return (int)discretization_.size();}
	int get_sizeDiscret_tilde()  const {return (int)discretization_tilde_.size();}
	Range get_range() const {return range_;}

	//! only for uniform discretization ...
	double get_discret_pas() const; 
	bool if_uniform() const;

	//! for extrapolation
	//void refine(int ratio);

	void print() const;
};

typedef boost::shared_ptr<Discretization_Simple> Discretization_Simple_PTR;
typedef boost::shared_ptr<const Discretization_Simple> Discretization_Simple_CONSTPTR;