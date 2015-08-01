#pragma once

#include <JBLMM/Instrument/GenericSwap.h> 

class GenericVanillaSwapPricer
{

public:
	//! constructor by default 

	// destructor
	virtual ~GenericVanillaSwapPricer(){}	

	//! To validate the result    initLibor[i] = L_i[T_0]
	double genericVanillaSwap_Analytical(GenericSwap_CONSTPTR genericVanillaSwap, const std::vector<double>& liborInitValue)const;
};

typedef boost::shared_ptr<GenericVanillaSwapPricer> GenericVanillaSwapPricer_PTR;
typedef boost::shared_ptr<const GenericVanillaSwapPricer> GenericVanillaSwapPricer_CONSTPTR;



