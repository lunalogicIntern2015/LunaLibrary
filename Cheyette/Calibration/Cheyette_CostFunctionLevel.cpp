#include "Cheyette_CostFunctionLevel.h"



//value: method to overload to compute the cost functon value in x.
//ici norme 2 = sqrt(sum of squares)
//m fixé, on fait varier sigma

//parametre array 1D
//returns difference between Swaption Market Quotes and Swaption Model Values	
QuantLib::Disposable<QuantLib::Array> Cheyette_CostFunctionLevel::values(const QuantLib::Array& param_sigma1D) const
{
	double marketQuote =  cheyetteMarketData_PTR_->getVolQuotes()[indexSwaption_] ;

	//setter le parametre sigma de cheyetteApproxPTR avec Array
	//std::cout << "DEBUG LEVEL : " << param_sigma1D[0] << std::endl ;
	cheyetteApprox_PTR_->updateLevel_calib(param_sigma1D[0], indexSwaption_ - 1) ; //maj de sigma + des buffers de l'approx
	
	VanillaSwaption_PTR pSwaption = cheyetteMarketData_PTR_->getVect_swaptions()[indexSwaption_] ;
	double modelQuote	= cheyetteApprox_PTR_->volBlack(pSwaption) ;	

	//std::cout << "vol Black : " << modelQuote << std::endl ;

	//array de dimension 1
	QuantLib::Array res(1) ;
	res[0] = modelQuote - marketQuote ;	

	return res ;
}

//	double modelQuote	= NumericalMethods::Black_impliedVolatility(modelPrice / annuity0, S0, strike, T) ;
//conversion prix -> vol
	//double T			= marketData_PTR_->get_aExpiry()[indexSwaption_] ;	//coTerminalSwaptionVol_PTR_->getVectorExpiry() ;
	//double strike		= marketData_PTR_->getStrikeATM()[indexSwaption_] ;	//coTerminalSwaptionVol_PTR_->getStrike() ;
	//double S0			= cheyetteApprox_PTR_->get_buffer_s0_() ;
	//double annuity0		= cheyetteApprox_PTR_->swapRateDenominator(0., 0., 0.) ;

	//double modelQuote	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
	//																		S0, strike, T) ;

	//std::cout << "volModel : " << modelQuote << ", volQuote : " << marketQuote << std::endl ;
	//o_ << "; volModel : ;" << modelQuote << "; volQuote : ;" << marketQuote << std::endl ;
