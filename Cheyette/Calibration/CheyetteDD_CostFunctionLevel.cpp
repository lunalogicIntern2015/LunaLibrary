#include "CheyetteDD_CostFunctionLevel.h"



//value: method to overload to compute the cost functon value in x.
//ici norme 2 = sqrt(sum of squares)
//m fixé, on fait varier sigma

//parametre array 1D
QuantLib::Disposable<QuantLib::Array> CheyetteDD_CostFunctionLevel::values(const QuantLib::Array& param_sigma1D) const
{
	// difference between Swaption Market Quotes and Swaption Model Values
	double marketQuote = coTerminalSwaptionVol_PTR_->getDiagonalSwaptionVol() ;
	
	//setter le parametre sigma de cheyetteApproxPTR avec Array
	cheyetteApprox_PTR_->updateSigma_calib(o_, param_sigma1D[0], indexSwaption_ - 1) ; //maj de sigma + des buffers de l'approx
	double modelPrice	= cheyetteApprox_PTR_->prixSwaptionApproxPiterbarg() ;
	//conversion prix -> vol
	double T			= coTerminalSwaptionVol_PTR_->getVectorExpiry() ;
	double strike		= coTerminalSwaptionVol_PTR_->getStrike() ;
	double S0			= cheyetteApprox_PTR_->get_buffer_s0_() ;
	double annuity0		= cheyetteApprox_PTR_->swapRateDenominator(0., 0., 0.) ;

	//double modelQuote	= NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike, T) ;

	double modelQuote	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
																			S0, strike, T) ;

	std::cout << "volModel : " << modelQuote << ", volQuote : " << marketQuote << std::endl ;
	o_ << "; volModel : ;" << modelQuote << "; volQuote : ;" << marketQuote << std::endl ;
	//array de dimension 1
	QuantLib::Array res(1) ;
	res[0] = modelQuote - marketQuote ;	

	return res ;
}
