#pragma once

#include <LMM/helper/LMMTenorStructure.h>

#include <LMM/LmmModel/Lmm.h>
#include <LMM/LmmModel/McLmm.h>   //Model
#include <LMM/LmmModel/McTerminalLmm.h>  

#include <JBLMM/Instrument/GeneticSwap.h>   //instrument
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/TargetCoupon.h>
#include <JBLMM/Pricer/McGeneticSwapLMMPricer.h>

class McGeneticTargetSwapLmmPricer	:public McGeneticSwapLMMPricer
{

public:
	//constructor et destructor
	McGeneticTargetSwapLmmPricer(const McLmm_PTR mcLmm)
		:McGeneticSwapLMMPricer(mcLmm)
	{};
	~McGeneticTargetSwapLmmPricer(){};

	virtual double swapNPV(GeneticSwap_CONSTPTR geneticSwap, size_t nbSimulation)const;

protected: 
	//! one simulation
	virtual double evaluateCouponLeg(	LMM::Index indexValuationDate,
										const CouponLeg_CONSTPTR couponLeg,
										const std::vector<double>& numeraire, 
										const matrix& liborMatrix,
										LMMTenorStructure_PTR lmmTenorStructure) const;
};
typedef boost::shared_ptr<McGeneticTargetSwapLmmPricer> McGeneticTargetSwapLmmPricer_PTR;
typedef boost::shared_ptr<const McGeneticTargetSwapLmmPricer> McGeneticTargetSwapLmmPricer_CONSTPTR;

