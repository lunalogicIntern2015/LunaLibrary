#pragma once

#include <Cheyette/Calibrator/CheyetteBaseCostFunction.h>

/*For a setup of your own function, you need to implement a class 
which derives from the CostFunction class and implements the following virtual functions :*/

class CheyetteDD_CostFunctionSkew : public CheyetteBaseCostFunction
{
public:
	CheyetteDD_CostFunctionSkew(	MarketData_PTR marketData_PTR, size_t indexSwaption,  
									CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR)
				: CheyetteBaseCostFunction(marketData_PTR, indexSwaption, cheyetteApprox_PTR)	{}

	virtual ~CheyetteDD_CostFunctionSkew()
	{ 
		//à compléter ?
	}

	//sigma fixé, on fait varier m pour prendre en compte le skew
	virtual Disposable<Array> values(const Array& param_m) const ;

	double volShift(size_t indexSwaption, double strike, double shift) const ;
};

typedef boost::shared_ptr<CheyetteDD_CostFunctionSkew> CheyetteDD_CostFunctionSkew_PTR;
typedef boost::shared_ptr<const CheyetteDD_CostFunctionSkew> CheyetteDD_CostFunctionSkew_CONSTPTR;

