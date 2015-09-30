#pragma once
#include <Cheyette/Pricer/MC_Cheyette.h>	
#include <Instrument/VanillaSwap.h>  
#include <Instrument/VanillaSwaption.h>

#include <LMM/helper/GenericPath.h>
#include <fstream>

class MC_Cheyette_VanillaSwaptionPricer : public MC_Cheyette
{

public:
	MC_Cheyette_VanillaSwaptionPricer(	CheyetteModel_PTR			cheyetteModel_PTR,
										RNGenerator_PTR				rnGenerator,
										LMMTenorStructure_PTR		pTenorStructure,
										double						fwdProbaT,
										size_t						discretizationBetweenDates   )
		:MC_Cheyette(cheyetteModel_PTR, rnGenerator, pTenorStructure, fwdProbaT, discretizationBetweenDates){}


	virtual ~MC_Cheyette_VanillaSwaptionPricer(){}

	//! Pricing at time T0=0
	//std::vector<double> : mean_x, IC_left, IC_right
	std::vector<double> price(VanillaSwaption_PTR pVanillaSwaption, size_t nbSimulation) const;
	

	std::vector<std::vector<double>> priceMultipleStrikes(VanillaSwaption_PTR pVanillaSwaption, 
															size_t nbSimulation, 
															std::vector<double> shifts_bp) const ;

	//pricing de swaption pour des strikes correspondant à une standardized moneyness dans [-5, 5]
		//res[0] = prixSwaptionsPlusieursStrikes ;
		//res[1] = IC_left ; 
		//res[2] = IC_right ; 
		//res[3] = strikes ; 
		//res[4] = moneyness ;
	std::vector<std::vector<double>> priceMultipleStrikes(	VanillaSwaption_PTR pVanillaSwaption, 
															size_t nbSimulation, 
															double S0, double sigma_ATM) const ;

	void print(VanillaSwaption_PTR vanillaSwaption, 
				std::vector<size_t> nbSimus, 
				std::vector<double> prixMC,
				std::vector<double> IC_inf,
				std::vector<double> IC_sup) const ;

	//void printMC_vs_approx(double approx, double b_barre, 
	//						double annuityA0, double swapRateS0, double volBlack, 
	//						double a, double b,
	//						VanillaSwaption_PTR vanillaSwaption, 
	//						std::vector<size_t> nbSimus, 
	//						std::vector<double> prixMC,
	//						std::vector<double> IC_inf,
	//						std::vector<double> IC_sup) const ;

	//void printMC_vs_approx(std::ofstream& o,
	//						double approx, double b_barre, 
	//						double annuityA0, double swapRateS0, double volBlack, 
	//						double a, double b,
	//						VanillaSwaption_PTR vanillaSwaption, 
	//						std::vector<size_t> nbSimus, 
	//						std::vector<double> prixMC,
	//						std::vector<double> IC_inf,
	//						std::vector<double> IC_sup) const ;
	
};


typedef boost::shared_ptr<MC_Cheyette_VanillaSwaptionPricer> MC_Cheyette_VanillaSwaptionPricer_PTR;
typedef boost::shared_ptr<const MC_Cheyette_VanillaSwaptionPricer> MC_Cheyette_VanillaSwaptionPricer_CONSTPTR;
