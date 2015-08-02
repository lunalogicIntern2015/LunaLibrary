#pragma once

#include "useful_function.h"

class PayoffYuan
{
public:
    virtual double get_payoff(double S) = 0;
	virtual double get_strike() = 0;
};

class Put: public PayoffYuan
{
	double K;
public:

	Put(double K_ip):K(K_ip){};

    double get_payoff(double S)
	{
	    return max(K-S,0);
	}

	double get_strike()
	{
		return K;
	}
};