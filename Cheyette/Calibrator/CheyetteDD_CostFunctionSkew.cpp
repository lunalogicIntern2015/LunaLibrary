#include "CheyetteDD_CostFunctionSkew.h"

//value: method to overload to compute the cost functon value in x.
//ici norme 2 = sqrt(sum of squares)
//m fixé, on fait varier sigma
Disposable<Array> CheyetteDD_CostFunctionSkew::values(const Array& param_m) const
{
	// difference between Swaption Market Quotes and Swaption Model Values
	std::vector<double> skewQuotes = coTerminalSwaptionSkew_PTR_->getDiagonalSwaptionSkew() ;
	size_t nbQuotes = skewQuotes.size()  ;  
	Array res(nbQuotes);

	for (size_t i = 0 ; i < nbQuotes ; ++i)
	{
		double skewQuote	= skewQuotes[i] ;
		//setter le parametre m de cheyetteApproxPTR avec Array
		cheyetteApprox_PTR_->updateM_calib(param_m) ; //maj de m + des buffers de l'approx

		double strikeATM	= coTerminalSwaptionSkew_PTR_->getStrike()[i] ;
		double shiftUp		= coTerminalSwaptionSkew_PTR_->getShift() ;
		double shiftDown	= - shiftUp ; 
		double volShiftUp	= volShift(i, strikeATM, shiftUp) ;
		double volShiftDown	= volShift(i, strikeATM, shiftDown) ;
		double skewModel	= (volShiftUp - volShiftDown) / (2 * shiftUp) ;
		std::cout << "skewModel : " << skewModel << ", skewQuote : " << skewQuote << std::endl ;
		res[i] = skewModel - skewQuote ;	
	}

	return res ;
}

//retourne la vol implicite pour ATM + shift bp
double CheyetteDD_CostFunctionSkew::volShift(size_t indexSwaption, double strike, double shift) const
{

		cheyetteApprox_PTR_->setStrike(strike + shift) ;
	
	//verification que strike est shifte
		std::cout	<< "strike apres shift : " << cheyetteApprox_PTR_->get_buffer_UnderlyingSwap_().get_strike() 
					<< ", vs : " << strike  +shift 
					<< std::endl ;
		double modelPrice	= cheyetteApprox_PTR_->prixSwaptionApproxPiterbarg() ;
	//conversion prix -> vol
		double T			= marketData_PTR_->get_aExpiry()[indexSwaption] ;
		double S0			= cheyetteApprox_PTR_->get_buffer_s0_() ;
		double modelQuote	= NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike + shift, T) ;
		
	//remise au strike (ATM)
		cheyetteApprox_PTR_->setStrike(strike) ;
	//verification que strike est shifte
		std::cout	<< "remise strike ATM : " << cheyetteApprox_PTR_->get_buffer_UnderlyingSwap_().get_strike() 
					<< ", vs : " << strike  
					<< std::endl ;
		std::cout << std::endl ;

		return modelQuote ;
}
