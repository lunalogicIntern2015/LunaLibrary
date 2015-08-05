#pragma once
#include <vector>
#include <boost/shared_ptr.hpp>

#include <LMM/helper/Name.h>
#include <Instrument/Instrument.h>

class CallableInstrument : public Instrument
{
protected:
	std::vector<LMM::Index> exerciseTimes_; 
public:
	//construtor
	CallableInstrument(const std::vector<LMM::Index>& exerciseTimes):exerciseTimes_(exerciseTimes){}

	//getter
	const std::vector<LMM::Index>& getExerciseTimes()const{return exerciseTimes_;}

	virtual Instrument_CONSTPTR getSubInstrument(const size_t indexStart, const size_t indexEnd) const = 0;
};
typedef boost::shared_ptr<CallableInstrument> CallableInstrument_PTR;
typedef boost::shared_ptr<const CallableInstrument> CallableInstrument_CONSTPTR;