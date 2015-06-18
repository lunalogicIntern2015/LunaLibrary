#pragma once

#include <vector>

#include <JBLMM/Instrument/GeneticSwap.h>
#include <JBLMM/Element/CouponLeg.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>

#include <LMM/helper/LMMTenorStructure.h>
#include <LMM/helper/TenorType.h>
#include <LMM/helper/Name.h>

class InstrumentFactory
{
public:
	static GeneticSwap_CONSTPTR createVanillaSwap(	double strike, 
													LMM::Index indexStart, 
													LMM::Index indexEnd, 
													Tenor floatingLegTenor, 
													Tenor fixedLegTenor,
													LMMTenorStructure_CONSTPTR swapStructure,
													double nominal);

	static GeneticSwap_CONSTPTR createStandardTARNSwap(	double strike,
														LMM::Index indexStart, 
												 		LMM::Index indexEnd, 
														Tenor floatingLegTenor, 
														Tenor fixedLegTenor,
														LMMTenorStructure_CONSTPTR swapStructure,
														double nominal,
														double target);

	static GeneticSwap_CONSTPTR createGeneticSwap(CouponLeg_CONSTPTR Leg1, CouponLeg_CONSTPTR Leg2);
};

