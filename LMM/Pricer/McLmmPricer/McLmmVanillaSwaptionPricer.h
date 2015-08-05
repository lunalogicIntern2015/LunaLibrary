
//! MClmm price Swap, don't need to be too general !!!  
#pragma once
#include <boost/shared_ptr.hpp>

#include <LMM/Mc/McLmm.h>        // Model
#include <Instrument/VanillaSwaption.h>  // Instrument
#include <LMM/Pricer/McLmmPricer/McLmmVanillaSwapPricer.h>



class McLmmVanillaSwaptionPricer : public McLmmVanillaSwapPricer // to reuse the already implemented code in VanillaSwapPricer.
{
public:
	McLmmVanillaSwaptionPricer(const McLmm_PTR& mcLmm)
		: McLmmVanillaSwapPricer(mcLmm)
	{};

	//! Pricing at time T0 = 0
	double price (const VanillaSwaption& vanillaSwaption, size_t nbSimulation)  const;
};

typedef boost::shared_ptr<McLmmVanillaSwaptionPricer> McLmmVanillaSwaptionPricer_PTR;