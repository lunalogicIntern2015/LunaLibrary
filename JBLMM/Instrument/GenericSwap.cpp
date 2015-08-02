#include "JBLMM/Instrument/GenericSwap.h"
#include <cstddef>

GenericSwap::GenericSwap(CouponLeg_CONSTPTR leg1, CouponLeg_CONSTPTR leg2)
	:
		leg1_(leg1),
		leg2_(leg2)
{
}

GenericSwap::GenericSwap(CouponLeg_CONSTPTR leg1, CouponLeg_CONSTPTR leg2, LMMTenorStructure_PTR lmmTenorStructure)
	:
		leg1_(leg1),
		leg2_(leg2),
		lmmTenorStructure_(lmmTenorStructure)
{
}

 Instrument_CONSTPTR GenericSwap::getSubGenericSwap(const size_t indexStart, const size_t indexEnd) const
{
	if(indexStart==indexEnd)
		return nullptr;
	return Instrument_CONSTPTR(new GenericSwap(	getLeg1()->getSubCouponLeg(indexStart,indexEnd), 
												getLeg2()->getSubCouponLeg(indexStart,indexEnd)));
}