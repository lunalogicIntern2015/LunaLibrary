#pragma once

#include <Cheyette/Calibration/CheyetteBaseCostFunction.h>


/*For a setup of your own function, you need to implement a class 
which derives from the CostFunction class and implements the following virtual functions :

*/

//calibration locale : indexSwaption
class Cheyette_CostFunctionLevel : public CheyetteBaseCostFunction
{
private:
	size_t indexSwaption_ ;		

public:

	Cheyette_CostFunctionLevel(	Cheyette_SwaptionPricer_Approx_PTR cheyetteApprox_PTR,
								CheyetteMarketData_PTR cheyetteMarketData_PTR,
								size_t indexSwaption)
		:	CheyetteBaseCostFunction(cheyetteApprox_PTR, cheyetteMarketData_PTR), indexSwaption_(indexSwaption){}

	virtual ~Cheyette_CostFunctionLevel(){}

	void setIndexSwaption(size_t indexSwaption){indexSwaption_ = indexSwaption ;}

	//m fixé, on fait varier sigma
	virtual QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param_sigma1D) const ;

	//print ici

};

typedef boost::shared_ptr<Cheyette_CostFunctionLevel>		Cheyette_CostFunctionLevel_PTR ;
typedef boost::shared_ptr<const Cheyette_CostFunctionLevel>	Cheyette_CostFunctionLevel_CONSTPTR;

