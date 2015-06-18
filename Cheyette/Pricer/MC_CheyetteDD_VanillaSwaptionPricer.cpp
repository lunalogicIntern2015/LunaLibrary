#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>



double MC_CheyetteDD_VanillaSwaptionPricer::price(	const VanillaSwaption& vanillaSwaption, 
													size_t nbSimulation)  const
{
	//std::vector<double> dates = mcCheyette_->getDatesOfSimulation() ;
	// 
	//VanillaSwap vanillaSwap = vanillaSwaption.getUnderlyingSwap() ;
	//std::vector<size_t> fixedLegIndexSchedule		= vanillaSwap.get_fixedLegPaymentIndexSchedule() ; 
	//std::vector<size_t> floatingLegIndexSchedule	= vanillaSwap.get_floatingLegPaymentIndexSchedule() ;
	//double fixed_tenor = vanillaSwap.get_fixedLegTenorType().YearFraction() ;
	//double float_tenor = vanillaSwap.get_floatingLegTenorType().YearFraction() ;
	//double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 
	//double strike = vanillaSwap.get_strike() ;
	//double fwdProbaT_ = mcCheyette_->getFwdProbaT() ;

	//double valueSwaption(0) ;
	//for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	//{
	//	mcCheyette_->simulate_Euler() ;
	//	std::vector<double> x_t_one_sim = mcCheyette_->get_x_t_Cheyette() ;
	//	std::vector<double> y_t_one_sim = mcCheyette_->get_y_t_Cheyette() ;

	//	//fixedLeg
	//	double fixedLegValue(0) ;
	//	for (size_t i = 0 ; i < fixedLegIndexSchedule.size() ; ++i)
	//	{
	//		double dateFlux = fixedLegIndexSchedule[i] * tenor_ref ;
	//		double flow = fixed_tenor * strike ; 
	//		size_t index = numeric::findClosestDate(dateFlux, dates) ;
	//		double x_t = x_t_one_sim[index] ;
	//		double y_t = y_t_one_sim[index] ;

	//		//double numeraireRatio = numeraire[indexValuationDate]/numeraire[paymentIndex];
	//		double numeraireRatio = mcCheyette_->getCheyetteDD_Model()->P(0, fwdProbaT_, 0, 0) / 
	//										mcCheyette_->getCheyetteDD_Model()->P(dateFlux, fwdProbaT_, x_t, y_t) ;

	//		fixedLegValue += flow * numeraireRatio ;
	//	}	
	//	//floatLeg
	//	double floatLegValue(0) ;
	//	for (size_t i = 0 ; i < floatingLegIndexSchedule.size() ; ++i)
	//	{
	//		//mettre variable avec fixingDat.year Fraction et paymentDate.yearFraction
	//		double dateFlux		= floatingLegIndexSchedule[i] * tenor_ref ;
	//		double dateFixing	= dateFlux - float_tenor ;
	//		size_t index = numeric::findClosestDate(dateFixing, dates) ;
	//		double x_t = x_t_one_sim[index] ;
	//		double y_t = y_t_one_sim[index] ;
	//		double flow = float_tenor * 
	//				mcCheyette_->getCheyetteDD_Model()->Libor(dateFixing,	dateFixing, dateFlux, x_t, y_t) ;

	//		//double numeraireRatio = numeraire[indexValuationDate]/numeraire[paymentIndex];
	//		double numeraireRatio = mcCheyette_->getCheyetteDD_Model()->P(0, fwdProbaT_, 0, 0) / 
	//										mcCheyette_->getCheyetteDD_Model()->P(dateFlux, fwdProbaT_, x_t, y_t) ;

	//		floatLegValue +=  flow * numeraireRatio ;
	//	}	

	//	valueSwaption += std::max(floatLegValue - fixedLegValue, 0.) ;
	//}
	////std::cout << "somme des flux positifs swaption : " << valueSwaption << std::endl;  
	//return valueSwaption / nbSimulation ;
	throw "non code" ;
	return 0 ;
}