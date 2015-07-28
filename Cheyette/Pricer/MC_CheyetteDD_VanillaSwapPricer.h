#pragma once
#include <Cheyette/Pricer/MC_Cheyette.h>	
#include <LMM/instrument/VanillaSwap.h>  


class MC_CheyetteDD_VanillaSwapPricer : public MC_Cheyette
{

public:
		MC_CheyetteDD_VanillaSwapPricer(	CheyetteDD_Model_PTR		cheyetteDD_Model,
											RNGenerator_PTR				rnGenerator,
											LMMTenorStructure_PTR		pTenorStructure,
											size_t						fwdProbaT,
											size_t						discretizationBetweenDates   )
		:MC_Cheyette(cheyetteDD_Model, rnGenerator, pTenorStructure, fwdProbaT, discretizationBetweenDates){}

	virtual ~MC_CheyetteDD_VanillaSwapPricer(){}

	//! Pricing at time T0=0
	//double swapRate(const VanillaSwap& vanillaSwap, size_t nbSimulation) const;
	double swapNPV (VanillaSwap_PTR vanillaSwap, size_t nbSimulation, 
					size_t valuationIndex,
					Tenor tenorFloatingLeg, Tenor tenorFixedLeg) const;


};


typedef boost::shared_ptr<MC_CheyetteDD_VanillaSwapPricer> MC_CheyetteDD_VanillaSwapPricer_PTR;
typedef boost::shared_ptr<const MC_CheyetteDD_VanillaSwapPricer> MC_CheyetteDD_VanillaSwapPricer_CONSTPTR;