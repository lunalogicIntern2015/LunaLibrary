#include "CheyetteBaseCostFunction.h"



QuantLib::Real CheyetteBaseCostFunction::value(const QuantLib::Array& param_array1D) const 
{
	QL_REQUIRE (param_array1D. size ()==1 , "Base Cost Function, param is 1-dim");

	QuantLib::Array diff_cost = values(param_array1D);	//retourne VolModel(paramArray) - VolMarket

	//Real res = 0;
	//for (size_t i = 0; i < diff_cost.size(); ++i)
	//{
	//	res += diff_cost[i]*diff_cost[i];
	//}

	//Real res = diff_cost[indexSwaption_] * diff_cost[indexSwaption_] ; //n'intervient qu'une seule swaption

	QuantLib::Real res = diff_cost[0] * diff_cost[0] ;
	res = sqrt(res);		
	return res;	
}