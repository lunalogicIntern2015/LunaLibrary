#pragma once
#include <Cheyette/Pricer/MC_Cheyette.h>	
#include <Instrument/VanillaSwap.h>  

class MC_Cheyette_VanillaSwapPricer : public MC_Cheyette
{

public:
		MC_Cheyette_VanillaSwapPricer(	CheyetteModel_PTR			cheyetteModel_PTR,
										RNGenerator_PTR				rnGenerator,
										LMMTenorStructure_PTR		pTenorStructure,
										double						fwdProbaT,
										size_t						discretizationBetweenDates)
		:MC_Cheyette(cheyetteModel_PTR, rnGenerator, pTenorStructure, fwdProbaT, discretizationBetweenDates){}

	virtual ~MC_Cheyette_VanillaSwapPricer(){}


	//! one simulation - pour les produits vanille derives
	double evaluateFloatLeg(	const size_t valuationIndex,
								const std::vector<size_t>& indexFloatLeg,
								const Tenor tenorFloatLeg) const;

	double evaluateFixedLeg(	const size_t valuationIndex,
								const std::vector<size_t>& indexFixedLeg,
								const Tenor tenorFixedLeg, 
								const double fixedRate) const;

	double swapNPV (const VanillaSwap_PTR pVanillaSwap, 
					const size_t nbSimulation, 
					const size_t valuationIndex = 0) const;

};


typedef boost::shared_ptr<MC_Cheyette_VanillaSwapPricer> MC_Cheyette_VanillaSwapPricer_PTR;
typedef boost::shared_ptr<const MC_Cheyette_VanillaSwapPricer> MC_Cheyette_VanillaSwapPricer_CONSTPTR;