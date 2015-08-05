#include <Instrument/GenericSwap/GenericSwaption.h>

GenericSwaption::GenericSwaption(const LMM::Index maturity, GenericSwap_CONSTPTR geneticSwap)
	:
	maturity_(maturity),
	genericSwap_(geneticSwap)
{
}

bool GenericSwaption::check()const		//suppose the start index is the smallest payment index
{
	LMM::Index leg1StartIndex=getGenericSwap()->getLeg1()->getLeg()[0]->getPaymentIndex();
	LMM::Index leg2StartIndex=getGenericSwap()->getLeg2()->getLeg()[0]->getPaymentIndex();
	return (getMaturity()<=leg1StartIndex&&getMaturity()<=leg2StartIndex);
}
