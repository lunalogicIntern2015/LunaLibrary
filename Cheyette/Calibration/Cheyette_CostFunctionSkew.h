#pragma once

#include <Cheyette/Calibration/CheyetteBaseCostFunction.h>


class Cheyette_CostFunctionSkew : public CheyetteBaseCostFunction
{
private:
	size_t indexSwaption_ ;	

public:
	Cheyette_CostFunctionSkew(	Cheyette_SwaptionPricer_Approx_PTR cheyetteApprox_PTR,
								CheyetteMarketData_PTR cheyetteMarketData_PTR,
								size_t indexSwaption)
		:	CheyetteBaseCostFunction(cheyetteApprox_PTR, cheyetteMarketData_PTR), indexSwaption_(indexSwaption){}

	virtual ~Cheyette_CostFunctionSkew(){}

	void setIndexSwaption(size_t indexSwaption){indexSwaption_ = indexSwaption ;}

	//les autres parametres étant fixés, on fait varier le paramètre de skew (m pour DD, b pour Quadratic)
	virtual QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param_m1D) const ;
	// ??
	// ??
	double volShift(double strike, double shift) const ;  //plus d'index swaption car il n'y a plus qu 'une swaption
};

typedef boost::shared_ptr<Cheyette_CostFunctionSkew> Cheyette_CostFunctionSkew_PTR;
typedef boost::shared_ptr<const Cheyette_CostFunctionSkew> Cheyette_CostFunctionSkew_CONSTPTR;

