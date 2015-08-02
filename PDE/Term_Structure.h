#pragma once
#include <iostream>
#include <math.h>

using namespace std;

class Term_Structure
{
public:
	virtual double get_ShortRate(double t)=0;
	virtual double get_DiscountFactor(double t)=0;
};


class Const_Term_Structure : public Term_Structure
{
public:
    double r;
    
	//! constructor! 
	Const_Term_Structure(double r_ip):r(r_ip){};

	double get_ShortRate(double t){return r;}
	double get_DiscountFactor(double t)
	{   
		//cout << "discount factor " << r << " , "<<  t << endl;
		return exp(-r*t);
	}
};