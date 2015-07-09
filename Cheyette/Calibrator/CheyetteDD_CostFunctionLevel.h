#pragma once

#include <Cheyette/Calibrator/CheyetteBaseCostFunction.h>

/*For a setup of your own function, you need to implement a class 
which derives from the CostFunction class and implements the following virtual functions :*/

class CheyetteDD_CostFunctionLevel : public CheyetteBaseCostFunction
{

public:
	CheyetteDD_CostFunctionLevel(	MarketData_PTR marketData_PTR, size_t indexSwaption,  
									CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR)
				: CheyetteBaseCostFunction(marketData_PTR, indexSwaption, cheyetteApprox_PTR)
	{}

	virtual ~CheyetteDD_CostFunctionLevel()
	{ 
		//à compléter ?
	}

	//m fixé, on fait varier sigma
	virtual Disposable<Array> values(const Array& param_sigma) const ;

};

typedef boost::shared_ptr<CheyetteDD_CostFunctionLevel> CheyetteDD_CostFunctionLevel_PTR;
typedef boost::shared_ptr<const CheyetteDD_CostFunctionLevel> CheyetteDD_CostFunctionLevel_CONSTPTR;

