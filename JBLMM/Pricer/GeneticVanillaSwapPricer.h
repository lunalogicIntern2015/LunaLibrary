#pragma once

#include <JBLMM/Instrument/GeneticSwap.h> 

class GeneticVanillaSwapPricer
{

public:
	//! constructor by default 

	// destructor
	virtual ~GeneticVanillaSwapPricer(){}	

	//! To validate the result    initLibor[i] = L_i[T_0]
	double geneticVanillaSwap_Analytical(GeneticSwap_CONSTPTR geneticVanillaSwap, const std::vector<double>& liborInitValue)const;
};

typedef boost::shared_ptr<GeneticVanillaSwapPricer> GeneticVanillaSwapPricer_PTR;
typedef boost::shared_ptr<const GeneticVanillaSwapPricer> GeneticVanillaSwapPricer_CONSTPTR;



