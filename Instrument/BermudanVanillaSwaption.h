#pragma once

#include <vector>

#include <LMM/Helper/Name.h>
#include <Instrument/VanillaSwap.h>

#include <Instrument/Instrument.h>
#include <Instrument/CallableOption/CallableInstrument.h>


class BermudanVanillaSwaption : public CallableInstrument
{
	VanillaSwap_CONSTPTR vanillaSwap_;

	//check: exerciseTikme in [fistIndex, lastindex]
	bool checkIfExerciseTime();
public:
	//getters
	VanillaSwap_CONSTPTR getVanillaSwap_PTR()const{return vanillaSwap_;}
	const std::vector<LMM::Index>& getExerciseTimes()const{return exerciseTimes_;}

	//constructor
	BermudanVanillaSwaption(VanillaSwap_CONSTPTR vanillaSwap,std::vector<LMM::Index>& exerciseTimes);

	//destructor
	virtual ~BermudanVanillaSwaption(){}

	//SubSwap's origin time is startIndex!!
	VanillaSwap_CONSTPTR getSubVanillaSwap(LMM::Index startIndex, LMM::Index endIndex)const;

	Instrument_CONSTPTR getSubInstrument(const LMM::Index startIndex, const LMM::Index endIndex)const;

	//check if exercise time is valuable
	//bool checkIfExerciseTimeValid()const;
};

typedef boost::shared_ptr<BermudanVanillaSwaption> BermudanVanillaSwaption_PTR;
typedef boost::shared_ptr<const BermudanVanillaSwaption> BermudanVanillaSwaption_CONSTPTR;

