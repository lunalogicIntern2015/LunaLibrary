#pragma once

#include <iostream>
#include <cassert>
#include <string>
#include <cassert>

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <ql/termstructures/volatility/abcd.hpp>

#include <LMM/helper/TypeConverter.h>
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/helper/LMMTenorStructure.h>

#include <LMM/LmmModel/Lmm.h>
#include <LMM/LmmModel/McLmm.h>   //Model
#include <LMM/LmmModel/McTerminalLmm.h>  

#include <JBLMM/Instrument/GenericSwap.h>   //instrument
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Pricer/McLmmPricer.h>


class McLmmGenericSwapPricer  : public McLmmPricer 
{

protected:

	McLmm_PTR mcLmm_; // model

public:

	McLmmGenericSwapPricer(const McLmm_PTR mcLmm)
		:mcLmm_(mcLmm)
	{};

	//gettor
	McLmm_PTR getMcLmm()const{return mcLmm_;}

	//! Pricing at time T0=0
	//double swapRate(const VanillaSwap& vanillaSwap, size_t nbSimulation) const;
	virtual double swapNPV(GenericSwap_CONSTPTR geneticSwap, size_t nbSimulation)const;

	//price for one simulation
	double price_on_oneSimalation(	Instrument_CONSTPTR instrument, 
									const matrix& liborMatrix, 
									const std::vector<double>& numeraire)const;
	//! 
	void resetGeneratorToinitSeed(){mcLmm_->get_RNGenerator()->resetGeneratorToinitSeed();}

protected: 

	//! one simulation
	virtual double evaluateCouponLeg(	 const LMM::Index indexValuationDate,
										 const CouponLeg_CONSTPTR couponLeg,
										 const std::vector<double>& numeraire, 
										 const matrix& liborMatrix,
										 LMMTenorStructure_CONSTPTR lmmTenorStructure) const;

	virtual double evaluateCappedFlooredCoupon(	CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
												double liborValue)const;

	double price_on_oneSimulation(	Instrument_CONSTPTR instrument, 
									LMM::Index indexValuationDate,
									const matrix& liborMatrix, 
									const std::vector<double>& numeraire) const;

};


typedef boost::shared_ptr<McLmmGenericSwapPricer> McLmmGenericSwapPricer_PTR;
typedef boost::shared_ptr<const McLmmGenericSwapPricer> McLmmGenericSwapPricer_CONSTPTR;
