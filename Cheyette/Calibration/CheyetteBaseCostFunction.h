#pragma once

#include <ql/types.hpp>
#include <ql/math/array.hpp>
#include <ql/math/optimization/costfunction.hpp>

#include <Cheyette/Calibration/CheyetteMarketData.h>

#include <Cheyette/Pricer/Cheyette_SwaptionPricer_Approx.h>
#include <Instrument/VanillaSwaption.h>


class CheyetteBaseCostFunction : public QuantLib::CostFunction
{
protected :
	mutable Cheyette_SwaptionPricer_Approx_PTR	cheyetteApprox_PTR_ ;
	CheyetteMarketData_PTR cheyetteMarketData_PTR_ ;

public:
	CheyetteBaseCostFunction(	Cheyette_SwaptionPricer_Approx_PTR cheyetteApprox_PTR,
								CheyetteMarketData_PTR cheyetteMarketData_PTR)
			: cheyetteApprox_PTR_(cheyetteApprox_PTR), cheyetteMarketData_PTR_(cheyetteMarketData_PTR) 
	{}

	virtual ~CheyetteBaseCostFunction(){}

	//value: method to overload to compute the cost functon value in x.
	//ici norme 2
	virtual QuantLib::Real value(const QuantLib::Array& param_array1D) const ;
	
	virtual QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param_sigma1D) const = 0 ;

	Cheyette_SwaptionPricer_Approx_PTR	getCheyette_ApproxPricer_PTR() const {return cheyetteApprox_PTR_ ;}
	CheyetteMarketData_PTR				getCheyetteMarketData_PTR() const {return cheyetteMarketData_PTR_ ;}



};

typedef boost::shared_ptr<CheyetteBaseCostFunction> CheyetteBaseCostFunction_PTR;
typedef boost::shared_ptr<const CheyetteBaseCostFunction> CheyetteBaseCostFunction_CONSTPTR;