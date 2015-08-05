#pragma once

#include <LMM/Helper/LMMTenorStructure.h>

#include <LMM/Model/Lmm.h>
#include <LMM/Mc/McLmm.h>   //Model
#include <LMM/Mc/McTerminalLmm.h>  

#include <Instrument/Instrument.h>
#include <Instrument/GenericSwap/GenericSwap.h>   //instrument
#include <Instrument/Coupon/Coupon.h>
#include <Instrument/Coupon/CappedFlooredCoupon.h>
#include <Instrument/Coupon/TargetCoupon.h>
#include <LMM/Pricer/McLmmPricer/McLmmGenericSwapPricer.h>

class McLmmGenericTargetSwapPricer	:public McLmmGenericSwapPricer
{

public:
	//constructor et destructor
	McLmmGenericTargetSwapPricer(const McLmm_PTR mcLmm)
		:McLmmGenericSwapPricer(mcLmm)
	{};
	//destructor
	~McLmmGenericTargetSwapPricer(){};

	double price_on_oneSimalation(	Instrument_CONSTPTR instrument, 
									LMM::Index indexValuationDate,
									const matrix& liborMatrix, 
									const std::vector<double>& numeraire) const;

	virtual double swapNPV(GenericSwap_CONSTPTR geneticSwap, size_t nbSimulation)const;

protected: 
	//! one simulation
	virtual double evaluateCouponLeg(	LMM::Index indexValuationDate,
										const CouponLeg_CONSTPTR couponLeg,
										const std::vector<double>& numeraire, 
										const matrix& liborMatrix,
										LMMTenorStructure_PTR lmmTenorStructure) const;
};
typedef boost::shared_ptr<McLmmGenericTargetSwapPricer> McLmmGenericTargetSwapPricer_PTR;
typedef boost::shared_ptr<const McLmmGenericTargetSwapPricer> McLmmGenericTargetSwapPricer_CONSTPTR;

