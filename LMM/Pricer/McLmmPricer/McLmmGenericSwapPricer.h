#pragma once

#include <iostream>
#include <cassert>
#include <string>
#include <cassert>

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <ql/termstructures/volatility/abcd.hpp>

#include <LMM/Helper/TypeConverter.h>
#include <RNGenerator/RNGenerator.h>
#include <LMM/Helper/LMMTenorStructure.h>

#include <LMM/Model/Lmm.h>
#include <LMM/Mc/McLmm.h>   //Model
#include <LMM/Mc/McTerminalLmm.h>  

#include <Instrument/GenericSwap/GenericSwap.h>   //instrument
#include <Instrument/Coupon/Coupon.h>
#include <Instrument/Coupon/CappedFlooredCoupon.h>
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>


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
