#include <cassert>
#include <iostream>
#include <cmath>

#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>



double MC_CheyetteDD_VanillaSwapPricer::swapNPV(double t_valo, 
												const VanillaSwap& vanillaSwap, 
												size_t nbSimulation)  const
{
	double dateFlux, x_t, y_t, fixedLegValue(0), floatLegValue(0), valueSwap(0) ;
	std::vector<double> x_t_one_sim, y_t_one_sim ;  //one simulation 
	std::vector<double> dates = mcCheyette_->getDatesOfSimulation_() ;

	std::vector<size_t> fixedLegIndexSchedule		= vanillaSwap.get_fixedLegPaymentIndexSchedule() ; 
	std::vector<size_t> floatingLegIndexSchedule	= vanillaSwap.get_floatingLegPaymentIndexSchedule() ;
	double fixed_tenor = vanillaSwap.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = vanillaSwap.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 
	double strike = vanillaSwap.get_strike() ;
	size_t index = numeric::findClosestDate(t_valo, dates) ;

	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		mcCheyette_->simulate_Euler() ;
		x_t_one_sim = mcCheyette_->get_x_t_Cheyette_() ;
		y_t_one_sim = mcCheyette_->get_y_t_Cheyette_() ;

			x_t = x_t_one_sim[index] ;
			y_t = y_t_one_sim[index] ;
		//fixedLeg
		for (size_t i = 0 ; i < fixedLegIndexSchedule.size() ; ++i)
		{
			dateFlux = fixedLegIndexSchedule[i] * tenor_ref ;
			fixedLegValue += fixed_tenor * strike * mcCheyette_->getCheyetteDD_Model_()->P(t_valo, dateFlux, x_t, y_t) ; 
		}	
		//floatLeg
		for (size_t i = 0 ; i < floatingLegIndexSchedule.size() ; ++i)
		{
			dateFlux = floatingLegIndexSchedule[i] * tenor_ref ;
			double libor = mcCheyette_->getCheyetteDD_Model_()->Libor(t_valo, dateFlux - float_tenor, dateFlux, x_t, y_t) ;
			floatLegValue += float_tenor * libor * mcCheyette_->getCheyetteDD_Model_()->P(t_valo, dateFlux, x_t, y_t) ; 
		}	
		valueSwap += floatLegValue - fixedLegValue ;
	}

	return valueSwap / nbSimulation ;
}


//double MC_CheyetteDD_VanillaSwapPricer:: swapRate(LMM::Index indexValuationDate,
//										 const VanillaSwap& vanillaSwap,
//										 const std::vector<double>& numeraire, 
//										 const std::vector<double>& xt_Cheyette,
//										 const std::vector<double>& yt_Cheyette) const
//{
//	double pvFloating = pvFloatingLeg( indexValuationDate,vanillaSwap,numeraire, xt_Cheyette, yt_Cheyette);
//	double annuity = pvFixedLeg   (indexValuationDate, vanillaSwap, numeraire);
//	return pvFloating / annuity;
//}