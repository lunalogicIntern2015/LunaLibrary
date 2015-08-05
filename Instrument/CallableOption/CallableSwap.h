#pragma once

#include <vector>

#include <boost/shared_ptr.hpp>

#include <LMM/helper/Name.h>

#include <JBLMM/Instrument/GenericSwap.h>
#include <JBLMM/Instrument/CallableInstrument.h>


class CallableGenericSwap : public CallableInstrument
{
	GenericSwap_CONSTPTR genericSwap_;
	//std::vector<LMM::Index> exerciseTimes_;
	//check: exerciseTikme in [fistIndex, lastindex]
	bool checkIfExerciseTime()const;
public:
	//getters
	GenericSwap_CONSTPTR getGenericSwap()const{return genericSwap_;}
	const std::vector<LMM::Index>& getExerciseTimes()const{return exerciseTimes_;}

	//constructor
	CallableGenericSwap(GenericSwap_CONSTPTR genetricSwap, const std::vector<LMM::Index>& exerciseTimes);

	//destructor
	virtual ~CallableGenericSwap(){}

	//check if exercise time is valuable

	Instrument_CONSTPTR getSubInstrument(const size_t indexStart, const size_t indexEnd) const
	{
		return getGenericSwap()->getSubGenericSwap(indexStart, indexEnd);
	}

	
};

typedef boost::shared_ptr<CallableGenericSwap> CallableGenericSwap_PTR;
typedef boost::shared_ptr<const CallableGenericSwap> CallableGenericSwap_CONSTPTR;
