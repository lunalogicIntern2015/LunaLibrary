#pragma once

#include <vector>
#include "Model.h"
#include <PDE/Discretization.h>
#include <string>
#include <PDE/Matrix/Matrix.h>
#include <boost/shared_ptr.hpp>


class PDE_2D_Model
{
public:
	virtual ~PDE_2D_Model(){};

	virtual double A_x(double t, int x_index, int y_index) const = 0;
	virtual double B_x(double t, int x_index, int y_index) const = 0;
	virtual double C_x(double t, int x_index, int y_index) const = 0;

	virtual double A_y(double t, int x_index, int y_index) const = 0;
	virtual double B_y(double t, int x_index, int y_index) const = 0;
	virtual double C_y(double t, int x_index, int y_index) const = 0;

	virtual double F_x_y(double t, int x_index, int y_index) const = 0;

	virtual std::string get_model_discret_type() const = 0;

	virtual void calculate_L(int t_index,  const Matrix& U_ip) const{};  // only for SLV
};

typedef boost::shared_ptr<PDE_2D_Model> PDE_2D_Model_PTR;
typedef boost::shared_ptr<const PDE_2D_Model> PDE_2D_Model_CONSTPTR;






