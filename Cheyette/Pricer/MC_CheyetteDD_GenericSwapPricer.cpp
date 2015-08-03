#include "MC_CheyetteDD_GenericSwapPricer.h"


//simulation
double MC_CheyetteDD_GenericSwapPricer::swapNPV(GenericSwap_CONSTPTR genericSwap, size_t nbSimulation, 
												Tenor tenorLeg1, Tenor tenorLeg2) const
{
	double TEST_sommeLeg1 = 0. ;
	double TEST_sommeLeg2 = 0. ;

	double result	= 0. ;
	size_t indexValuationDate = 0 ;

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimulation == 0){std::cout << double(itrSimulation)/double(nbSimulation)*100 << "%" << std::endl ;}
		//simulate_Euler() ;
		//double npv1  = evaluateCouponLeg(indexValuationDate, genericSwap->getLeg1(), tenorLeg1);

		//double npv2  = evaluateCouponLeg(indexValuationDate, genericSwap->getLeg2(), tenorLeg2);			
		//TEST_sommeLeg1 += npv1 ;
		//TEST_sommeLeg2 += npv2 ;
		//result += npv1 - npv2;
	}
	double TEST_meanLeg1 = TEST_sommeLeg1 / nbSimulation ;
	double TEST_meanLeg2 = TEST_sommeLeg2 / nbSimulation ;
	result   /=nbSimulation; 

	std::cout   << "prix MC swap : " << result << std::endl;
	std::cout   << "prix leg1 : " << TEST_meanLeg1 << std::endl;
	std::cout   << "prix leg2 : " << TEST_meanLeg2 << std::endl;
	
	return result;
}


