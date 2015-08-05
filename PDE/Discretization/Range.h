#pragma once

#include <iostream>
#include <boost/shared_ptr.hpp>
//! D=1
class Range
{
private:
	//! interval [left,right], require leftBorder_ < rightBorder_
	double leftBorder_; 
	double rightBorder_; 

public:
	//! constructor
	Range(double x1, double x2);
	virtual ~Range(){};
	double get_leftBorder(){return leftBorder_;}
	double get_rightBorder(){return rightBorder_;}
};


typedef boost::shared_ptr<Range> Range_PTR;
typedef boost::shared_ptr<const Range> Range_CONSTPTR;