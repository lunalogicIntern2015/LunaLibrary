#include <cassert>
#include <iostream>
#include <cmath>

#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>


//simulation
double MC_CheyetteDD_VanillaSwapPricer::swapNPV (VanillaSwap_PTR vanillaSwap, size_t nbSimulation, 
												 size_t valuationIndex,
												 Tenor tenorFloatingLeg, Tenor tenorFixedLeg) const
{
	assert(tenorFloatingLeg.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(tenorFixedLeg.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(vanillaSwap->get_EndDate() <= fwdProbaT_) ;

	double TEST_sommeLeg1 = 0. ;
	double TEST_sommeLeg2 = 0. ;

	double result	= 0. ;
	double result_2	= 0. ;

	std::vector<size_t> floatingLegIndexSchedule	= vanillaSwap->get_floatingLegPaymentIndexSchedule() ;
	std::vector<size_t> fixedLegIndexSchedule		= vanillaSwap->get_fixedLegPaymentIndexSchedule() ; 

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimulation == 0){std::cout << double(itrSimulation)/double(nbSimulation)*100 << "%" << std::endl ;}
		simulate_Euler() ;
		double npv1  = evaluateFloatLeg(valuationIndex, floatingLegIndexSchedule, tenorFloatingLeg);

		double npv2  = evaluateFixedLeg(valuationIndex, fixedLegIndexSchedule, tenorFixedLeg, vanillaSwap->get_strike());

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

