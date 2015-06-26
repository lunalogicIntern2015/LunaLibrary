#include <cassert>
#include <iostream>
#include <cmath>

#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>


//simulation
double MC_CheyetteDD_VanillaSwapPricer::swapNPV (VanillaSwap_PTR vanillaSwap, size_t nbSimulation, 
												 Tenor tenorFloatingLeg, Tenor tenorFixedLeg) const
{
	double TEST_sommeLeg1 = 0. ;
	double TEST_sommeLeg2 = 0. ;

	double result	= 0. ;
	size_t indexValuationDate = 0 ;

	std::vector<size_t> floatingLegIndexSchedule	= vanillaSwap->get_floatingLegPaymentIndexSchedule() ;
	std::vector<size_t> fixedLegIndexSchedule		= vanillaSwap->get_fixedLegPaymentIndexSchedule() ; 

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		simulate_Euler() ;
		double npv1  = evaluateFloatLeg(indexValuationDate, floatingLegIndexSchedule, numeraires_, 
										x_t_Cheyette_, y_t_Cheyette_, tenorFloatingLeg);

		double npv2  = evaluateFixedLeg(indexValuationDate, fixedLegIndexSchedule, numeraires_, 
										x_t_Cheyette_, y_t_Cheyette_, tenorFixedLeg, vanillaSwap->get_strike());

		TEST_sommeLeg1 += npv1 ;
		TEST_sommeLeg2 += npv2 ;
		result += npv1 - npv2;
	}
	double TEST_meanLeg1 = TEST_sommeLeg1 / nbSimulation ;
	double TEST_meanLeg2 = TEST_sommeLeg2 / nbSimulation ;
	result   /=nbSimulation; 

	std::cout   << "prix MC swap : " << result << std::endl;
	std::cout   << "prix leg1 : " << TEST_meanLeg1 << std::endl;
	std::cout   << "prix leg2 : " << TEST_meanLeg2 << std::endl;
	
	return result;
}

