#pragma once

#include <LMM/helper/LMMTenorStructure.h>

#include <LMM/LmmModel/Lmm.h>
#include <LMM/LmmModel/McLmm.h>   //Model
#include <LMM/LmmModel/McTerminalLmm.h>  

#include <JBLMM/Instrument/Instrument.h>
#include <JBLMM/Instrument/GenericSwap.h>   //instrument
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/TargetCoupon.h>
#include <JBLMM/Pricer/McLmmGenericSwapPricer.h>

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

