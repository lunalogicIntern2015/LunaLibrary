#include "Cheyette_CostFunctionConvexity.h"

void Cheyette_CostFunctionConvexity::calculateConvexity()
{
	std::vector<double> vect_vol_convexity_pp = cheyetteMarketData_2_PTR_->getVect_vol_convexity_pp() ;
	size_t nbQuotesITM = vect_vol_convexity_pp.size() ;
	std::vector<double> vect_vol_convexity_mm = cheyetteMarketData_2_PTR_->getVect_vol_convexity_mm() ;
	size_t nbQuotesOTM = vect_vol_convexity_mm.size() ;


	std::vector<double> vect_vol_ATM = cheyetteMarketData_PTR_->getVolQuotes() ; //format : 0.2 pour 20% -- pas des %, ni des bp. 
	size_t nbQuotesATM = vect_vol_ATM.size() ;				

	assert(nbQuotesITM == nbQuotesOTM) ;
	assert(nbQuotesITM == nbQuotesATM) ;

	//quotes en %, shift en %
	double shift = cheyetteMarketData_2_PTR_->getShiftConvexity() / 100 ;
	
	for (size_t i = 0 ; i < nbQuotesITM ; ++i)
	{
		double vectVolATM_pourcent = vect_vol_ATM[i] * 100. ;
		double convexity = (vect_vol_convexity_mm[i]/100. 
							- 2 * vectVolATM_pourcent /100.
							+ vect_vol_convexity_pp[i] /100.) / (shift * shift) ;
		convexity_.push_back(convexity) ;
	}
	convexity_.shrink_to_fit() ;
}

//value: method to overload to compute the cost functon value in x.
//calculs sur données en %
QuantLib::Disposable<QuantLib::Array> Cheyette_CostFunctionConvexity::values(const QuantLib::Array& param_c_1D) const
{

	double convexityQuote = convexity_[indexSwaption_] ; 
	//std::cout << "DEBUG : convexity mkt value : " << convexityQuote << std::endl ;
	//std::cout << "DEBUG : param convexity  : " << param_c_1D << std::endl ;
	cheyetteApprox_PTR_->updateConvexity_calib(param_c_1D[0], indexSwaption_ - 1) ; 

	double strikeATM	= cheyetteMarketData_PTR_->getStrikeATM()[indexSwaption_] /* * 100 */ ;	
	VanillaSwaption_PTR pSwaption = cheyetteMarketData_PTR_->getVect_swaptions()[indexSwaption_] ;
	double volATM		= cheyetteApprox_PTR_->volBlack(pSwaption) /* * 100 */ ;

	double shiftUp		= cheyetteMarketData_2_PTR_->getShiftConvexity() / 100. ;	
	double shiftDown	= - shiftUp ; 

	assert(strikeATM > shiftUp) ;  //K_ATM -> K_ATM +/- shift

	double volShiftUp	= volShift(strikeATM, shiftUp) /* * 100 */ ;  
	//std::cout << "DEBUG : vol shift up : " << volShiftUp << std::endl ;
	double volShiftDown	= volShift(strikeATM, shiftDown) /* * 100 */ ;
	//std::cout << "DEBUG : vol shift down : " << volShiftDown << std::endl ;

	double convexityModel	= (volShiftDown - 2. * volATM + volShiftUp) / (shiftUp * shiftUp) ;
	//std::cout << "DEBUG : volATM : " << volATM << std::endl ;
	//std::cout << "DEBUG : shiftUp : " << shiftUp << std::endl ;
	//std::cout << "DEBUG : convexityModel : " << convexityModel << std::endl ;

	QuantLib::Array res(1) ;
	res[0] = convexityModel - convexityQuote ;		
	return res ;
}

//retourne la vol implicite pour ATM + shift bp
double Cheyette_CostFunctionConvexity::volShift(double strike, double shift) const
{
	VanillaSwaption_PTR pSwaption = cheyetteMarketData_PTR_->getVect_swaptions()[indexSwaption_] ;
	pSwaption->set_strike(strike + shift) ;

	double modelQuote	= cheyetteApprox_PTR_->volBlack(pSwaption) ;
	
	//remise au strike (ATM)
	pSwaption->set_strike(strike) ;
		
	return modelQuote ;
}


