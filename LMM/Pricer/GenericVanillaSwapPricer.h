#pragma once

#include <Instrument/GenericSwap/GenericSwap.h> 

//! only for the test purpose, need to delete it latter ! 
//! this is to validate the genericVanillaSwap instrument

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



