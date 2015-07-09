#pragma once

#include <deque>

#include <ql/types.hpp>
#include <ql/math/array.hpp>
#include <ql/math/optimization/costfunction.hpp>

#include <Cheyette/Calibrator/MarketData.h>
#include <Cheyette/Calibrator/CoTerminalSwaptionQuotes.h>

#include <Cheyette/Pricer/CheyetteDD_VanillaSwaptionApproxPricer.h>
#include <LMM/instrument/VanillaSwaption.h>

using namespace QuantLib ;  //pour CostFunction, Real, Array


class CheyetteBaseCostFunction : public CostFunction
{
protected :
// 	CoTerminalSwaptionQuotes_CONSTPTR					coTerminalSwaptionQuotes_PTR_ ;

	mutable size_t										indexSwaption_ ; //sur quel index de swaption on calibre 								
	mutable CheyetteDD_VanillaSwaptionApproxPricer_PTR	cheyetteApprox_PTR_ ;
	
public:
	CheyetteBaseCostFunction(	size_t indexSwaption,  
								CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR)
				:	indexSwaption_(indexSwaption), 
					cheyetteApprox_PTR_(cheyetteApprox_PTR)
	{}

	virtual ~CheyetteBaseCostFunction()
	{ 
		//à compléter ?
	}

	//value: method to overload to compute the cost functon value in x.
	//ici norme 2 = sqrt(sum of squares	
	virtual Real value(const Array& param_array) const ;
	
	virtual Disposable<Array> values(const Array& param_sigma) const = 0 ;

	size_t										getIndexSwaption() const {return indexSwaption_ ;}
	CheyetteDD_VanillaSwaptionApproxPricer_PTR	getCheyetteDD_ApproxPricer_PTR() const {return cheyetteApprox_PTR_ ;}

	//void print(const std::string& filename) const;

};

typedef boost::shared_ptr<CheyetteBaseCostFunction> CheyetteBaseCostFunction_PTR;
typedef boost::shared_ptr<const CheyetteBaseCostFunction> CheyetteBaseCostFunction_CONSTPTR;