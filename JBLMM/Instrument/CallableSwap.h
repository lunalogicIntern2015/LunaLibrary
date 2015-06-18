#pragma once

#include <vector>

#include <JBLMM/Instrument/GeneticSwap.h>
#include <LMM/helper/Name.h>

class CallableSwap
{
	GeneticSwap_CONSTPTR geneticSwap_;
	std::vector<LMM::Index> exerciseTimes_;
	//check: exerciseTikme in [fistIndex, lastindex]
public:
	//getters
	GeneticSwap_CONSTPTR getGeneticSwap()const{return geneticSwap_;}
	const std::vector<LMM::Index>& getExerciseTimes()const{return exerciseTimes_;}

	//constructor
	CallableSwap(GeneticSwap_CONSTPTR geneticSwap, const std::vector<LMM::Index>& exerciseTimes);

	//destructor
	virtual ~CallableSwap(){}

	//check if exercise time is valuable
	//bool checkIfExerciseTimeValid()const;
	
};

