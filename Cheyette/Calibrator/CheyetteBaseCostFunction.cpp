#include "CheyetteBaseCostFunction.h"



Real CheyetteBaseCostFunction::value(const Array& param_array) const 
{
	Array diff_cost = values(param_array);	//retourne VolModel(paramArray) - VolMarket

	//Real res = 0;
	//for (size_t i = 0; i < diff_cost.size(); ++i)
	//{
	//	res += diff_cost[i]*diff_cost[i];
	//}

	Real res = diff_cost[indexSwaption_] * diff_cost[indexSwaption_] ; //à vérifier, n'intervient qu'une seule swaption

	res = sqrt(res);		
	return res;	
}