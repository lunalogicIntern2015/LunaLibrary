#pragma once
#include <Cheyette/Pricer/MC_Cheyette.h>	
#include <Instrument/GenericSwap/GenericSwaption.h>

class MC_Cheyette_GenericSwaptionPricer : public MC_Cheyette
{

public:
	MC_Cheyette_GenericSwaptionPricer(	CheyetteModel_PTR			cheyetteModel_PTR,
										RNGenerator_PTR				rnGenerator,
										LMMTenorStructure_PTR		pTenorStructure,
										size_t						fwdProbaT,
										size_t						discretizationBetweenDates   )
		:MC_Cheyette(cheyetteModel_PTR, rnGenerator, pTenorStructure, fwdProbaT, discretizationBetweenDates){}

	virtual ~MC_Cheyette_GenericSwaptionPricer(){}

	//! Pricing at time T0=0
	std::vector<double> price(GenericSwaption_CONSTPTR genericSwaption, size_t nbSimulation) const;

	void print(GenericSwaption_CONSTPTR genericSwaption, 
				std::vector<size_t> nbSimus, 
				std::vector<double> prixMC,
				std::vector<double> IC_inf,
				std::vector<double> IC_sup) const ;

	void printMC_vs_approx(double approx, double b_barre, 
							double annuityA0, double swapRateS0, double volBlack, 
							double a, double b,
							GenericSwaption_CONSTPTR genericSwaption, 
							std::vector<size_t> nbSimus, 
							std::vector<double> prixMC,
							std::vector<double> IC_inf,
							std::vector<double> IC_sup) const ;



};

typedef boost::shared_ptr<MC_Cheyette_GenericSwaptionPricer>       MC_Cheyette_GenericSwaptionPricer_PTR;
typedef boost::shared_ptr<const MC_Cheyette_GenericSwaptionPricer> MC_Cheyette_GenericSwaptionPricer_CONSTPTR;