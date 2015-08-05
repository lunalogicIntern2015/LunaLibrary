#include <Instrument/BermudanVanillaSwaption.h>
#include <LMM/Helper/LMMTenorStructure.h>
#include <cstddef>

BermudanVanillaSwaption::BermudanVanillaSwaption(VanillaSwap_CONSTPTR vanillaSwap, std::vector<LMM::Index>& exerciseTimes)
	:
	CallableInstrument(exerciseTimes),
	vanillaSwap_(vanillaSwap)
{	
	if(!checkIfExerciseTime())
		throw("CallableSwap constructor, checkIfExerciseTime() failed.");
}

VanillaSwap_CONSTPTR BermudanVanillaSwaption::getSubVanillaSwap(LMM::Index startIndex, LMM::Index endIndex)const
{
	VanillaSwap_CONSTPTR vanillaswap_PTR = getVanillaSwap_PTR();
	assert(startIndex<endIndex);
	assert(startIndex>=vanillaswap_PTR->get_indexStart()&&endIndex<=vanillaswap_PTR->get_indexEnd());
	double strike = vanillaswap_PTR->get_strike();

	//copy pointer or deep copy for Tenor or LMMTenorStruture? 
	const Tenor subFloatingTenor(vanillaswap_PTR->get_floatingLegTenorType());
	const Tenor subFixedTenor(vanillaswap_PTR->get_fixedLegTenorType());

	size_t floatingFixedRatio = subFixedTenor.ratioTo(subFloatingTenor);
	assert((endIndex-startIndex)%floatingFixedRatio==0);

	LMMTenorStructure_PTR LTS = vanillaswap_PTR->get_LMMTenorStructure();
	LMMTenorStructure_PTR subLTS(new LMMTenorStructure(*LTS.get()));

	return VanillaSwap_CONSTPTR(new VanillaSwap(	strike,
													startIndex,
													endIndex,
													subFloatingTenor,
													subFixedTenor,
													subLTS
												));
}

bool BermudanVanillaSwaption::checkIfExerciseTime()
{
	std::sort(exerciseTimes_.begin(),exerciseTimes_.end());
	LMM::Index startIndex = vanillaSwap_->get_indexStart();
	LMM::Index endIndex = vanillaSwap_->get_indexEnd();
	if(exerciseTimes_[0]<startIndex && exerciseTimes_.back()>endIndex)
		return false;

	return true;
}

Instrument_CONSTPTR BermudanVanillaSwaption::getSubInstrument(const LMM::Index startIndex, const LMM::Index endIndex)const
{

	if(startIndex==endIndex)
		return nullptr;

	VanillaSwap_CONSTPTR vanillaswap_PTR = getVanillaSwap_PTR();
	assert(startIndex<endIndex);
	assert(startIndex>=vanillaswap_PTR->get_indexStart()&&endIndex<=vanillaswap_PTR->get_indexEnd());
	double strike = vanillaswap_PTR->get_strike();

	//copy pointer or deep copy for Tenor or LMMTenorStruture? 
	const Tenor subFloatingTenor(vanillaswap_PTR->get_floatingLegTenorType());
	const Tenor subFixedTenor(vanillaswap_PTR->get_fixedLegTenorType());

	size_t floatingFixedRatio = subFixedTenor.ratioTo(subFloatingTenor);
	assert((endIndex-startIndex)%floatingFixedRatio==0);

	LMMTenorStructure_PTR LTS = vanillaswap_PTR->get_LMMTenorStructure();
	LMMTenorStructure_PTR subLTS(new LMMTenorStructure(*LTS.get()));

	return Instrument_CONSTPTR(new VanillaSwap(	strike,
													startIndex,
													endIndex,
													subFloatingTenor,
													subFixedTenor,
													subLTS
												));	
}


