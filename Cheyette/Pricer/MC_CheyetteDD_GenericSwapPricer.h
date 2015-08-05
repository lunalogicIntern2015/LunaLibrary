#pragma once
#include <Cheyette/Pricer/MC_Cheyette.h>	
#include <Instrument/GenericSwap/GenericSwap.h>

class MC_CheyetteDD_GenericSwapPricer : public MC_Cheyette
{

public:
	MC_CheyetteDD_GenericSwapPricer(	CheyetteDD_Model_PTR		cheyetteDD_Model,
										RNGenerator_PTR				rnGenerator,
										LMMTenorStructure_PTR		pTenorStructure,
										size_t						fwdProbaT,
										size_t						discretizationBetweenDates   )
		:MC_Cheyette(cheyetteDD_Model, rnGenerator, pTenorStructure, fwdProbaT, discretizationBetweenDates){}

	virtual ~MC_CheyetteDD_GenericSwapPricer(){}


	//! Pricing at time T0=0
	double swapNPV(GenericSwap_CONSTPTR genericSwap, size_t nbSimulation, Tenor tenorLeg1, Tenor tenorLeg2) const ;



};

typedef boost::shared_ptr<MC_CheyetteDD_GenericSwapPricer>       MC_CheyetteDD_GenericSwapPricer_PTR;
typedef boost::shared_ptr<const MC_CheyetteDD_GenericSwapPricer> MC_CheyetteDD_GenericSwapPricer_CONSTPTR;