#include <cassert>
#include <iostream>
#include <cmath>

#include <Cheyette/Pricer/MC_Cheyette_VanillaSwapPricer.h>

//! one simulation - pour les produits vanille 
double MC_Cheyette_VanillaSwapPricer::evaluateFloatLeg(	const size_t valuationIndex,
														const std::vector<size_t>& indexFloatLeg,
														const Tenor tenorFloatLeg) const
{
	//verifie que grille TenorStructure plus fine que les flux du swap
	double tenorFloatYearFrac	= tenorFloatLeg.YearFraction() ;
	double tenorStructYearFrac	= pTenorStructure_->get_tenorType().YearFraction() ;
	assert(tenorStructYearFrac <= tenorFloatYearFrac) ;

	double price = 0.0;

	//rapport tenorFloat / tenorStructure  (ex : swap 6M sur TenorStructure 3M) rapport de 2 = nbIndices
	size_t nbIndices = static_cast<size_t>(tenorFloatYearFrac / tenorStructYearFrac) ;

	for(size_t i=0; i<indexFloatLeg.size(); ++i)
	{
		size_t paymentIndex		= indexFloatLeg[i] ;

		if(paymentIndex <= valuationIndex)continue;

		double paymentDate		= paymentIndex * tenorStructYearFrac ;
		size_t fixingIndex		= paymentIndex - nbIndices ;
		double fixingDate		= fixingIndex * tenorStructYearFrac ;	

		//  !!! FAUX !!!
		//double libor = cheyetteModel_PTR_->libor(fixingDate, fixingDate, paymentDate, 
		//										x_t_Cheyette_[fixingIndex], y_t_Cheyette_[fixingIndex]) ; 

		double libor = cheyetteModel_PTR_->libor(valuationIndex * tenorStructYearFrac, fixingDate, paymentDate, 
												x_t_Cheyette_[valuationIndex], y_t_Cheyette_[valuationIndex]) ; 	

		//payoffFlow = nominal * deltaFloat * value  (ici NOMINAL = 1)
		double payoffFlow		= tenorFloatYearFrac * libor ;		//tenorFloat = deltaFloat																			
		double numeraireRatio	= numeraires_[valuationIndex]/numeraires_[paymentIndex];
		price += payoffFlow * numeraireRatio;
	}

//version 2 : mieux mais encore un petit ecart

// A MULTIPLIER PAR LES NUMERAIRES !!!
	//double startDateSwap	= (indexFloatLeg[0] - nbIndices ) * tenorStructYearFrac ;
	//double endDateSwap		= indexFloatLeg[indexFloatLeg.size() - 1] * tenorStructYearFrac ;
	//double price = 1 - cheyetteDD_Model_->P(startDateSwap, endDateSwap, x_t_Cheyette_[valuationIndex], 
	//																	y_t_Cheyette_[valuationIndex]) ;

	return price;
}


double MC_Cheyette_VanillaSwapPricer::evaluateFixedLeg(	const size_t valuationIndex,
										const std::vector<size_t>& indexFixedLeg,
										const Tenor tenorFixedLeg, 
										const double fixedRate) const
{
	//verifie que grille TenorStructure plus fine que les flux du swap
	double tenorFixedYearFrac	= tenorFixedLeg.YearFraction() ;
	double tenorStructYearFrac	= pTenorStructure_->get_tenorType().YearFraction() ;
	assert(tenorStructYearFrac <= tenorFixedYearFrac) ;

	double price = 0.0;
	for(size_t i=0; i<indexFixedLeg.size(); ++i)
	{
		size_t paymentIndex		= indexFixedLeg[i] ;

		if(paymentIndex <= valuationIndex)continue;

		double paymentDate		= paymentIndex * tenorStructYearFrac ;
			
		//payoffFlow = nominal * deltaFloat * value  (ici NOMINAL = 1)
		double payoffFlow		= tenorFixedYearFrac * fixedRate ;
																					
		double numeraireRatio	= numeraires_[valuationIndex]  / numeraires_[paymentIndex];
		price += payoffFlow*numeraireRatio;
	}

	return price;
}


double MC_Cheyette_VanillaSwapPricer::swapNPV (	const VanillaSwap_PTR pVanillaSwap, 
												const size_t nbSimulation, 
												const size_t valuationIndex) const
{
	Tenor tenorFloatingLeg =  pVanillaSwap->get_floatingLegTenorType() ;
	Tenor tenorFixedLeg = pVanillaSwap->get_fixedLegTenorType() ;

	assert(tenorFloatingLeg.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(tenorFixedLeg.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(pVanillaSwap->get_EndDate() <= fwdProbaT_) ;

	double TEST_sommeLeg1 = 0. ;
	double TEST_sommeLeg2 = 0. ;

	double result	= 0. ;
	double result_2	= 0. ;

	std::vector<size_t> floatingLegIndexSchedule	= pVanillaSwap->get_floatingLegPaymentIndexSchedule() ;
	std::vector<size_t> fixedLegIndexSchedule		= pVanillaSwap->get_fixedLegPaymentIndexSchedule() ; 

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimulation == 0){std::cout << double(itrSimulation)/double(nbSimulation)*100 << "%" << std::endl ;}
		simulate_Euler() ;
		double npv1  = evaluateFloatLeg(valuationIndex, floatingLegIndexSchedule, tenorFloatingLeg);

		double npv2  = evaluateFixedLeg(valuationIndex, fixedLegIndexSchedule, tenorFixedLeg, pVanillaSwap->get_strike());

		TEST_sommeLeg1 += npv1 ;
		TEST_sommeLeg2 += npv2 ;
		double res = npv1 - npv2;
		result += res ;
		result_2 += res* res ;
	}
	double TEST_meanLeg1 = TEST_sommeLeg1 / nbSimulation ;
	double TEST_meanLeg2 = TEST_sommeLeg2 / nbSimulation ;

	double mean_x	= result / nbSimulation; 
	double mean_x2	= result_2 / nbSimulation; 
 
	double variance = mean_x2 - mean_x * mean_x ;

	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimulation);
	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimulation);

	std::cout   << "prix MC swap : " << mean_x << std::endl;
	std::cout   << "prix leg1 : " << TEST_meanLeg1 << std::endl;
	std::cout   << "prix leg2 : " << TEST_meanLeg2 << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;

	return mean_x ;
}

