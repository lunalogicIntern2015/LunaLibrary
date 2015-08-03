#pragma once

#include <deque>

#include <ql/types.hpp>
#include <ql/math/array.hpp>
#include <ql/math/optimization/costfunction.hpp>

#include <Cheyette/Calibrator/MarketData.h>
#include <Cheyette/Calibrator/CoTerminalSwaptionQuotes.h>

#include <Cheyette/Pricer/CheyetteDD_VanillaSwaptionApproxPricer.h>
#include <LMM/instrument/VanillaSwaption.h>

//using namespace QuantLib ;  //pour CostFunction, Real, Array


class CheyetteBaseCostFunction : public QuantLib::CostFunction
{
protected :
	size_t indexSwaption_ ;		
	mutable CheyetteDD_VanillaSwaptionApproxPricer_PTR	cheyetteApprox_PTR_ ;
	ostream& o_ ;
public:
	CheyetteBaseCostFunction(	ostream& o, CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR, 
								size_t indexSwaption)
			: o_(o), cheyetteApprox_PTR_(cheyetteApprox_PTR), indexSwaption_(indexSwaption)
	{}

	virtual ~CheyetteBaseCostFunction(){}

	//value: method to overload to compute the cost functon value in x.
	//ici norme 2 = sqrt(sum of squares	
	virtual QuantLib::Real value(const QuantLib::Array& param_array1D) const ;
	
	virtual QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param_sigma1D) const = 0 ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR	getCheyetteDD_ApproxPricer_PTR() const {return cheyetteApprox_PTR_ ;}

	//void print(const std::string& filename) const;

};

typedef boost::shared_ptr<CheyetteBaseCostFunction> CheyetteBaseCostFunction_PTR;
typedef boost::shared_ptr<const CheyetteBaseCostFunction> CheyetteBaseCostFunction_CONSTPTR;