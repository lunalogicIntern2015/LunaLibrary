#pragma once

#include <vector>

#include <LMM/helper/LMMTenorStructure.h>
#include <LMM/helper/TenorType.h>
#include <LMM/helper/Name.h>

#include <JBLMM/Instrument/GenericSwap.h>
#include <JBLMM/Element/CouponLeg.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>
#include <JBLMM/Element/ConstRate.h>



class InstrumentFactory
{
public:
	static GenericSwap_CONSTPTR createVanillaSwap(	double strike, 
													LMM::Index indexStart, 
													LMM::Index indexEnd, 
													Tenor floatingLegTenor, 
													Tenor fixedLegTenor,
													LMMTenorStructure_CONSTPTR swapStructure,
													double nominal);

	static GenericSwap_CONSTPTR createStandardTARNSwap(	double strike,
														LMM::Index indexStart, 
												 		LMM::Index indexEnd, 
														Tenor floatingLegTenor, 
														Tenor fixedLegTenor,
														LMMTenorStructure_CONSTPTR swapStructure,
														double nominal,
														double target);

	static GenericSwap_CONSTPTR createGenericSwap(CouponLeg_CONSTPTR Leg1, CouponLeg_CONSTPTR Leg2);
};

