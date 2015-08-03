//#pragma once
//
//#include <iostream>
//#include "Range.h"
//#include <vector>
//
//
//class Discretization_D2
//{
//	//std::vector<double> discretization_t_;
//	//std::vector<double> discretization_x_;
//	//std::vector<double> discretization_y_;
//
//	//! this is only for the uniform discret
//	//double Discret_t_pas_;   // = range/sizeDiscretization
//	//double Discret_x_pas_;
//	//double Discret_y_pas_;
//
//public:
//
//	virtual std::vector<double> get_discret_t() = 0;
//	virtual std::vector<double> get_discret_x() = 0;
//	virtual std::vector<double> get_discret_y() = 0;
//	virtual std::vector<double> get_discret_x_tilde() = 0;
//	virtual std::vector<double> get_discret_y_tilde() = 0;
//
//	virtual double get_discret_t	    (int i_index) const = 0;
//	virtual double get_discret_x		(int i_index) const = 0;
//	virtual double get_discret_x_tilde  (int i_index) const = 0;
//	virtual double get_discret_y        (int i_index) const = 0;
//	virtual double get_discret_y_tilde  (int i_index) const = 0;
//
//	virtual int get_sizeDiscret_t() = 0; 
//	virtual int get_sizeDiscret_x() = 0; 
//	virtual int get_sizeDiscret_y() = 0; 
//	virtual int get_sizeDiscret_x_tilde() = 0; 
//	virtual int get_sizeDiscret_y_tilde() = 0; 
//
//
//	//! this is only for uniform_discretization
//	//double get_sizeDiscret_t_pas(){return Discret_t_pas_;}
//	//double get_sizeDiscret_x_pas(){return Discret_x_pas_;}
//	//double get_sizeDiscret_y_pas(){return Discret_y_pas_;}
//
//	virtual double get_discret_t(int i) = 0;
//	virtual double get_discret_x(int i) = 0;
//	virtual double get_discret_y(int i) = 0;
//	virtual double get_discret_x_tilde(int i) = 0;
//	virtual double get_discret_y_tilde(int i) = 0;
//
//
//	virtual double get_discret_t_left_border()  = 0; 
//	virtual double get_discret_t_right_border() = 0;
//	virtual double get_discret_x_left_border(double t)  = 0; 
//	virtual double get_discret_x_right_border(double t) = 0; 
//	virtual double get_discret_y_left_border(double t)  = 0; 
//	virtual double get_discret_y_right_border(double t) = 0; 
//
//};
//
