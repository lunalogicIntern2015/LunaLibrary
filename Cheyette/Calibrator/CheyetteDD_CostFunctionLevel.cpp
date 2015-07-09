#include "CheyetteDD_CostFunctionLevel.h"



//value: method to overload to compute the cost functon value in x.
//ici norme 2 = sqrt(sum of squares)
//m fixé, on fait varier sigma
Disposable<Array> CheyetteDD_CostFunctionLevel::values(const Array& param_sigma) const
{
	// difference between Swaption Market Quotes and Swaption Model Values
	size_t nbQuotes = marketData_PTR_->get_aExpiry().size() ;
	Array res(nbQuotes);
	std::vector<double> volQuotes = marketData_PTR_->getVolQuotes() ;

	for (size_t i = 0 ; i < nbQuotes ; ++i)
	{
		double marketQuote	= volQuotes[i] ;
		//setter le parametre sigma de cheyetteApproxPTR avec Array
		cheyetteApprox_PTR_->updateSigma_calib(param_sigma) ; //maj de sigma + des buffers de l'approx
		double modelPrice	= cheyetteApprox_PTR_->prixSwaptionApproxPiterbarg() ;
		//conversion prix -> vol
		double T			= marketData_PTR_->get_aExpiry()[i] ;
		double strike		= marketData_PTR_->getStrikeATM()[i] ;
		double S0			= cheyetteApprox_PTR_->get_buffer_s0_() ;
		double modelQuote	= NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike, T) ;
	
		res[i] = modelQuote - marketQuote ;	
	}

	return res ;
}
