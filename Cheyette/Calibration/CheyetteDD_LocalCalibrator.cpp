#include "CheyetteDD_LocalCalibrator.h"

#include <fstream>
#include <cassert>


//minimisations 1D : sigma, puis m, autant de fois que nbIterations_
//calibration sur UNE swaption
void CheyetteDD_LocalCalibrator::solve() const
{
	assert(nbIterations_ >= 1) ;
	QuantLib::LevenbergMarquardt	minimizationSolver ; //(functionEpsilon(), rootEpsilon(), 1e-16);

	//minimisation sur sigma
	std::vector<double> vSigma1D(1) ; vSigma1D[0] = calibrated_sigma_[currentSwaptionIndex_] ;
	QuantLib::Array sigma1D = HelperArray::vectorToArray(vSigma1D) ;
	minimizePositiveConstraint(sigma1D, minimizationSolver, costFunctionLevel_PTR_) ;  

	for (size_t i = 0 ; i < nbIterations_ ; ++i)
	{
		//minimisation sur m
		std::vector<double> vM1D(1) ; vM1D[0] = calibrated_m_[currentSwaptionIndex_] ;
		QuantLib::Array m1D = HelperArray::vectorToArray(vM1D) ;
		minimizeBoundaryConstraint(m1D, minimizationSolver, costFunctionSkew_PTR_) ; 
		
		//minimisation sur sigma
		std::vector<double> vSigma1D(1) ; vSigma1D[0] = calibrated_sigma_[currentSwaptionIndex_] ;
		QuantLib::Array sigma1D = HelperArray::vectorToArray(vSigma1D) ;
		minimizePositiveConstraint(sigma1D, minimizationSolver, costFunctionLevel_PTR_) ;
	}

	// !!! A VOIR 
	//mise à jour de calibrated_m_sigma et calibrated_m ok car passage par référence dans minimizeConstraint

////pour le smile

//	double strike		= costFunctionLevel_PTR_->getCoTerminalSwaptionVol_PTR()->getStrike() ;
//	double shift		= costFunctionSkew_PTR_->getCoTerminalSwaptionSkew_PTR()->getShift() ;
//
//	o_ << " ; vol ATM - 5bp ; vol ATM ; vol ATM + 5bp " << std::endl ;
//	o_ << "STRIKE ; " << strike - shift << " ; " << strike << " ; " << strike + shift << std::endl ;
////vol ATM - shift modele 
//	double T			= costFunctionLevel_PTR_->getCoTerminalSwaptionVol_PTR()->getVectorExpiry() ;
//	double S0			= costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->get_buffer_s0_() ;
//	double annuity0		= costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->swapRateDenominator(0., 0., 0.) ;
//
//	costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->setStrike(strike - shift) ;
//	double modelPrice	= costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->prixSwaptionApprox() ;
//	double volATMmoins	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
//																			S0, strike - shift, T) ;
//	o_ << "VOL MODEL ; " << volATMmoins ;
//
////vol ATM 
//	costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->setStrike(strike) ;
//	modelPrice	= costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->prixSwaptionApprox() ;
//	double volATM	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
//																			S0, strike, T) ;
//	o_ << " ; " << volATM ;
//
////vol ATM + shift modele 
//	costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->setStrike(strike + shift) ;
//	modelPrice	= costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->prixSwaptionApprox() ;
//	double volATMplus	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
//																			S0, strike + shift, T) ;
//	o_ << " ; " << volATMplus << std::endl ;
//
//
//	costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->setStrike(strike) ;
//
	//costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->get_CheyetteDD_Model()->show() ;
}

//boucle sur toutes les swaptions de calibration
void CheyetteDD_LocalCalibrator::calibrate() const
{
	CheyetteMarketData_PTR marketData = costFunctionSkew_PTR_->getCheyetteMarketData_PTR() ;
	std::vector<VanillaSwaption_PTR> pVectSwaptions = marketData->getVect_swaptions() ;
	//nb de swaptions coterminales
	size_t nbSwaptionsCalibration = pVectSwaptions.size() - 1 ; //swaption à l'indice 0 ne compte pas

	assert(marketData == costFunctionSkew_PTR_->getCheyetteMarketData_PTR() ) ;
	assert(nbSwaptionsCalibration > 0) ;

	for (size_t indexSwaption = 1 ; indexSwaption <= nbSwaptionsCalibration ; ++indexSwaption)
	{
		std::cout << endl ;
		std::cout << "------ Calibration sur la swaption " ;
		//			<< indexSwaption << "Y" << coterminal - indexSwaption << "Y -------" << std::endl ;
		std::cout << std::endl ;

		currentSwaptionIndex_ = indexSwaption ; 

		//initialisation de Cheyette APPROX sur la swaption en cours
		costFunctionLevel_PTR_->getCheyette_ApproxPricer_PTR()->preCalculateALL(pVectSwaptions[currentSwaptionIndex_]) ;

		costFunctionLevel_PTR_->setIndexSwaption(indexSwaption) ;
		costFunctionSkew_PTR_->setIndexSwaption(indexSwaption) ;

		solve() ;

		//o << "VOL MARKET ; "
		//	<< volATMmm_ptr->get_UpperTriangularVanillaSwaptionQuotes()(indexSwaption, coterminal - indexSwaption).second	<< ";"
		//	<< volATM_ptr->get_UpperTriangularVanillaSwaptionQuotes()(indexSwaption, coterminal - indexSwaption).second		<< ";"
		//	<< volATMpp_ptr->get_UpperTriangularVanillaSwaptionQuotes()(indexSwaption, coterminal - indexSwaption).second	<<std::endl ;
		//o << endl ;
	}
	costFunctionSkew_PTR_->getCheyette_ApproxPricer_PTR()->get_CheyetteModel()->show() ;  //print(o) ;
}