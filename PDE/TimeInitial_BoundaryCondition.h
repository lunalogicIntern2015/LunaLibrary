#pragma once

#include "Term_Structure.h"
#include "Payoff.h"
#include <boost/shared_ptr.hpp>

class TimeInitial_BoundaryCondition
{
public:
	 virtual ~TimeInitial_BoundaryCondition(){};

	 virtual double get_inital_bc(double x, double y) const= 0;  
};
typedef boost::shared_ptr<TimeInitial_BoundaryCondition>       TimeInitial_BoundaryCondition_PTR;
typedef boost::shared_ptr<const TimeInitial_BoundaryCondition> TimeInitial_BoundaryCondition_CONSTPTR;


