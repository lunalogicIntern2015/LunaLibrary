#pragma once

#include <vector>

#include <LMM/Helper/LMMTenorStructure.h>
#include <LMM/Helper/TenorType.h>
#include <LMM/Helper/Name.h>

#include <Instrument/GenericSwap/GenericSwap.h>
#include <Instrument/Coupon/CouponLeg.h>
#include <Instrument/Coupon/CappedFlooredCoupon.h>
#include <Instrument/Rate/LiborRate.h>
#include <Instrument/Rate/ConstRate.h>



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

