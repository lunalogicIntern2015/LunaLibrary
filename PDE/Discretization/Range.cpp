#include <PDE/Discretization/Range.h>

Range::Range(double x1, double x2)
	{
		if(x1 == x2)
		{
			std::cout << "error of x1==x2, in Range::Range(double x1, double x2) " << std::endl; // throw exeption later 
		}
		else if(x1<x2)
		{
		    leftBorder_ = x1;
            rightBorder_ = x2;
		}
		else
		{
		    leftBorder_ = x2;
            rightBorder_ = x1;
		}
	}
