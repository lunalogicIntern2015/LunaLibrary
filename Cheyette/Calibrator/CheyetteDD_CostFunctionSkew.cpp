#include "CheyetteDD_CostFunctionSkew.h"

//value: method to overload to compute the cost functon value in x.
//ici norme 2 = sqrt(sum of squares)
Disposable<Array> CheyetteDD_CostFunctionSkew::values(const Array& param_m1D) const
{
	// difference between Swaption Market Quotes and Swaption Model Values
	double skewQuote = coTerminalSwaptionSkew_PTR_->getDiagonalSwaptionSkew() ;

	//setter le parametre m de cheyetteApproxPTR avec Array
	cheyetteApprox_PTR_->updateM_calib(o_, param_m1D[0], indexSwaption_ - 1) ; //maj de m + des buffers de l'approx

	double strikeATM	= coTerminalSwaptionSkew_PTR_->getStrike() ;
	double shiftUp		= coTerminalSwaptionSkew_PTR_->getShift() ;
	double shiftDown	= - shiftUp ; 
	double volShiftUp	= volShift(strikeATM, shiftUp) ;
	double volShiftDown	= volShift(strikeATM, shiftDown) ;
	double skewModel	= (volShiftUp - volShiftDown) / (2 * shiftUp) ;
	std::cout << "skewModel : " << skewModel << ", skewQuote : " << skewQuote << std::endl ;
	o_ << "; skewModel : ;" << skewModel << "; skewQuote : ;" << skewQuote << std::endl ;

	Array res(1) ;
	res[0] = skewModel - skewQuote ;		
	return res ;
}

//retourne la vol implicite pour ATM + shift bp
double CheyetteDD_CostFunctionSkew::volShift(double strike, double shift) const
{
	cheyetteApprox_PTR_->setStrike(strike + shift) ;	
	//verification que strike est shifte
		//std::cout	<< "strike apres shift : " << cheyetteApprox_PTR_->get_buffer_UnderlyingSwap_().get_strike() 
		//			<< ", vs : " << strike  +shift 
		//			<< std::endl ;
	double modelPrice	= cheyetteApprox_PTR_->prixSwaptionApproxPiterbarg() ;
	//conversion prix -> vol
	double T			= coTerminalSwaptionSkew_PTR_->getVectorExpiry() ;
	double S0			= cheyetteApprox_PTR_->get_buffer_s0_() ;
	double annuity0		= cheyetteApprox_PTR_->swapRateDenominator(0., 0., 0.) ;

	//double modelQuote	= NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike + shift, T) ;
	
	double modelQuote	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
																			S0, strike + shift, T) ;

	//remise au strike (ATM)
	cheyetteApprox_PTR_->setStrike(strike) ;
	//verification que strike est shifte
		//std::cout	<< "remise strike ATM : " << cheyetteApprox_PTR_->get_buffer_UnderlyingSwap_().get_strike() 
		//			<< ", vs : " << strike  
		//			<< std::endl ;
		//std::cout << std::endl ;

	return modelQuote ;
}
