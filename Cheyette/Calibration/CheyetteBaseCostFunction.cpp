#include "CheyetteBaseCostFunction.h"



QuantLib::Real CheyetteBaseCostFunction::value(const QuantLib::Array& param_array1D) const 
{
	QL_REQUIRE (param_array1D. size ()==1 , "Base Cost Function, param is 1-dim");

	QuantLib::Array diff_cost = values(param_array1D);	//retourne VolModel(paramArray) - VolMarket

	QuantLib::Real res = diff_cost[0] * diff_cost[0] ;
	res = sqrt(res);		
	return res;	
}

