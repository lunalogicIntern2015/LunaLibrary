#pragma once
#include <Cheyette/Pricer/MC_Cheyette.h>	
#include <LMM/instrument/VanillaSwap.h>  
#include <LMM/instrument/VanillaSwaption.h>


class MC_CheyetteDD_VanillaSwaptionPricer : public MC_Cheyette
{

public:
	MC_CheyetteDD_VanillaSwaptionPricer(CheyetteDD_Model_PTR		cheyetteDD_Model,
										RNGenerator_PTR				rnGenerator,
										LMMTenorStructure_PTR		pTenorStructure,
										size_t						fwdProbaT,
										size_t						discretizationBetweenDates   )
		:MC_Cheyette(cheyetteDD_Model, rnGenerator, pTenorStructure, fwdProbaT, discretizationBetweenDates){}


	virtual ~MC_CheyetteDD_VanillaSwaptionPricer(){}

	//! Pricing at time T0=0
	std::vector<double> price(VanillaSwaption_PTR pVanillaSwaption, size_t nbSimulation) const;

	void print(VanillaSwaption_PTR vanillaSwaption, 
				std::vector<size_t> nbSimus, 
				std::vector<double> prixMC,
				std::vector<double> IC_inf,
				std::vector<double> IC_sup) const ;

	void printMC_vs_approx(double approx, double b_barre, 
							double annuityA0, double swapRateS0, double volBlack, 
							double a, double b,
							VanillaSwaption_PTR vanillaSwaption, 
							std::vector<size_t> nbSimus, 
							std::vector<double> prixMC,
							std::vector<double> IC_inf,
							std::vector<double> IC_sup) const ;

	void printMC_vs_approx(std::ofstream& o,
							double approx, double b_barre, 
							double annuityA0, double swapRateS0, double volBlack, 
							double a, double b,
							VanillaSwaption_PTR vanillaSwaption, 
							std::vector<size_t> nbSimus, 
							std::vector<double> prixMC,
							std::vector<double> IC_inf,
							std::vector<double> IC_sup) const ;
	
};


typedef boost::shared_ptr<MC_CheyetteDD_VanillaSwaptionPricer> MC_CheyetteDD_VanillaSwaptionPricer_PTR;
typedef boost::shared_ptr<const MC_CheyetteDD_VanillaSwaptionPricer> MC_CheyetteDD_VanillaSwaptionPricer_CONSTPTR;
