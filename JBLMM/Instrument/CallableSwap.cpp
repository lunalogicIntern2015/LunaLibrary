#include "JBLMM/Instrument/CallableSwap.h"


CallableGenericSwap::CallableGenericSwap(GenericSwap_CONSTPTR genericSwap, const std::vector<LMM::Index>& exerciseTimes)
	:
	CallableInstrument(exerciseTimes),
	genericSwap_(genericSwap)
{
	if(!checkIfExerciseTime())
		throw("CallableSwap constructor, checkIfExerciseTime() failed.");
}

bool CallableGenericSwap::checkIfExerciseTime()const
{
	bool res = true;
	for(size_t i=0; i<exerciseTimes_.size(); i++)
	{
		LMM::Index startIndex	=	std::min(	genericSwap_->getLeg1()->getLeg()[0]->getPaymentIndex(), 
												genericSwap_->getLeg2()->getLeg()[0]->getPaymentIndex());
		LMM::Index endIndex		=	std::max(	genericSwap_->getLeg1()->getLeg().back()->getPaymentIndex(), 
												genericSwap_->getLeg2()->getLeg().back()->getPaymentIndex());
		if(exerciseTimes_[i] <	startIndex-1 || exerciseTimes_[i]>endIndex)
		{
			res=false;
			break;
		}
	}
	return res;		
}