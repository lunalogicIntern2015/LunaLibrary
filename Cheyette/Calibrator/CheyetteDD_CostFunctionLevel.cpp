#include "CheyetteDD_CostFunctionLevel.h"



//value: method to overload to compute the cost functon value in x.
//ici norme 2 = sqrt(sum of squares)
//m fixé, on fait varier sigma
Disposable<Array> CheyetteDD_CostFunctionLevel::values(const Array& param_sigma) const
{
	// difference between Swaption Market Quotes and Swaption Model Values
	std::vector<double> volQuotes = coTerminalSwaptionVol_PTR_->getDiagonalSwaptionVol() ;
	size_t nbQuotes = volQuotes.size()  ;  
	Array res(nbQuotes);

	for (size_t i = 0 ; i < nbQuotes ; ++i)
	{
		double marketQuote	= volQuotes[i] ;
		//setter le parametre sigma de cheyetteApproxPTR avec Array
		cheyetteApprox_PTR_->updateSigma_calib(param_sigma) ; //maj de sigma + des buffers de l'approx
		double modelPrice	= cheyetteApprox_PTR_->prixSwaptionApproxPiterbarg() ;
		//conversion prix -> vol
		double T			= coTerminalSwaptionVol_PTR_->getVectorExpiry()[i] ;
		double strike		= coTerminalSwaptionVol_PTR_->getStrike()[i] ;
		double S0			= cheyetteApprox_PTR_->get_buffer_s0_() ;
		double modelQuote	= NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike, T) ;
	
		res[i] = modelQuote - marketQuote ;	
	}

	return res ;
}
