#include "Cheyette_CostFunctionSkew.h"


QuantLib::Disposable<QuantLib::Array> Cheyette_CostFunctionSkew::values(const QuantLib::Array& param_m1D) const
{
	// difference between Swaption Market Quotes and Swaption Model Values
	double skewQuote = cheyetteMarketData_PTR_->getSkew()[indexSwaption_] ;  

	//setter le parametre m de cheyetteApproxPTR avec Array

	//std::cout << "DEBUG SKEW : " << param_m1D[0] << std::endl ;
	cheyetteApprox_PTR_->updateSkew_calib(param_m1D[0], indexSwaption_ - 1) ; 

	double strikeATM	= cheyetteMarketData_PTR_->getStrikeATM()[indexSwaption_] ;		
	double shiftUp		= cheyetteMarketData_PTR_->getShift() ;							
	double shiftDown	= - shiftUp ; 
	double volShiftUp	= volShift(strikeATM, shiftUp) ;
	double volShiftDown	= volShift(strikeATM, shiftDown) ;
	double skewModel	= (volShiftUp - volShiftDown) / (2 * shiftUp) ;

	//std::cout << "skewModel : " << skewModel << ", skewQuote : " << skewQuote << std::endl ;

	QuantLib::Array res(1) ;
	res[0] = skewModel - skewQuote ;		
	return res ;
}

//retourne la vol implicite pour ATM + shift bp
double Cheyette_CostFunctionSkew::volShift(double strike, double shift) const
{
	VanillaSwaption_PTR pSwaption = cheyetteMarketData_PTR_->getVect_swaptions()[indexSwaption_] ;
	pSwaption->set_strike(strike + shift) ;

	//verification que strike est shifte
		//std::cout	<< "strike apres shift : " << cheyetteApprox_PTR_->get_buffer_UnderlyingSwap_().get_strike() 
		//			<< ", vs : " << strike  +shift << std::endl ;

	double modelQuote	= cheyetteApprox_PTR_->volBlack(pSwaption) ;
	
	//remise au strike (ATM)
	pSwaption->set_strike(strike) ;
		//std::cout	<< "remise strike ATM : " << cheyetteApprox_PTR_->get_buffer_UnderlyingSwap_().get_strike() 
		//			<< ", vs : " << strike  << std::endl ;
		
	return modelQuote ;
}
