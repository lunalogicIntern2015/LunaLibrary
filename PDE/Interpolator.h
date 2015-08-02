#pragma once 

#include <iostream>
#include <string>
#include <PDE/Matrix/Matrix.h>

using namespace std;

class Interpolator       //!  1.D 
{
public: 
	vector<double> x_;  //! x_grid
	vector<double> y_;  //! y_value
	vector<double> y2;  //! second derivative of y, for intermediate calculation.
	int n_;
	int jsav_;
	int dj_;

	//! default constructor do nothing 
	Interpolator();
	Interpolator(const vector<double>& x_grid, const vector<double>& y_value);

	virtual ~Interpolator(){};

	double  interpolate(double x);
	double  calculate_integral_approximation() const;

	int     locate(double x);
	void    set_y2();
	vector<double> interpolate(vector<double>& x);
};


//class Interpolator_D2
//{
//public:
//	vector<double> x_grid_;   //! x_grid size:s1
//	vector<double> y_grid_;   //! y_grid size:s2
//    vector<vector<double>> grid_value_;  //! size:s1*s2
//
//public:
//	Interpolator_D2(const vector<double>& x_grid, const vector<double>& y_grid, const vector<vector<double>>& value);
//	double  interpolate(double x, double y);
//	double  calculate_integral_approximation() const;
//	
//	int     locate(double x, double y);
//	void    set_2nd_deirvative_x();
//	void    set_2nd_deirvative_y();
//
//	vector<double> interpolate(vector<vector<double>>& grid);
//};


class Interpolator_D2
{
public:
	vector<double> x_grid_;   //! x_grid size:s1
	vector<double> y_grid_;   //! y_grid size:s2
    Matrix grid_value_;  //! size:s1*s2

	//Interpolator		 x_interp;
	vector<Interpolator> y_interp_v;

public:

	Interpolator_D2(const vector<double>& x_grid,
					const vector<double>& y_grid,
					const Matrix& value);

	double  interpolate(double x, double y){throw("Error Interpolator_D2::interpolate(double x, double y) not implemented!");}
	Matrix interpolate(const vector<double>& x_grid, const vector<double>& y_grid);
	//double  calculate_integral_approximation() const;
	
	//int     locate(double x, double y);
	//void    set_2nd_deirvative_x();
	//void    set_2nd_deirvative_y();

	//vector<double> interpolate(vector<vector<double>>& grid);
};

