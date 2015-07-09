#include "CheyetteDD_LocalCalibrator.h"

#include <fstream>
#include <cassert>
#include <string.h>
#include <ctime>


//void CheyetteLocalCalibrator::retrieve_calib_info()
//{
//	retrieve_calib_global_error();
//
//	boost::shared_ptr<LmmLocal_gCostFunction> pLmmLocal_gCostFunction = boost::static_pointer_cast<LmmLocal_gCostFunction>(pLmmBaseCostFunction_);
//
//	calibrated_gDelegate_matrix_ = pLmmLocal_gCostFunction->get_GMatrixMapping()->get_gDelegate_Ref();
//
//	const size_t swaption_matrix_size = calibrated_gDelegate_matrix_.size1();
//
//	size_t nb_errors_counter=0;
//	for(size_t iExperity=1;iExperity<swaption_matrix_size;++iExperity)
//	{
//		const size_t maxTenorBound =  swaption_matrix_size - iExperity;
//
//		for(size_t jTenor=1;jTenor<maxTenorBound;++jTenor)
//		{
//			// bad calibrated param
//			if(calibrated_gDelegate_matrix_(iExperity,jTenor) <1e-4
//				//&& calibrated_gDelegate_matrix_(iExperity,jTenor) > -5000 // case truncated, values are setted  
//					)
//			{
//				std::pair<size_t,size_t> bad_param_indices(iExperity,jTenor);
//				bad_swaptions_indices_.push_back(bad_param_indices);
//			}
//
//			// Compute param error
//			if(isVirtualCalibration_)
//			{				
//				const double relative_param_error = 100*(calibrated_gDelegate_matrix_(iExperity,jTenor) - true_gDelegate_matrix_(iExperity,jTenor))/true_gDelegate_matrix_(iExperity,jTenor);
//				rel_error_gDelegate_matrix_(iExperity,jTenor) = relative_param_error;
//				if(std::abs(relative_param_error) > max_param_rel_error_ ) max_param_rel_error_ = std::abs(relative_param_error);
//			}
//		}
//	}
//
//	// ---------------- PRINT GENERAL RESULT INFO OF CALBRATION  ---------------------------
//	std::stringstream general_info;
//
//	Shifted_HGVolatilityParam::ABCDParameter abcdParam = pLmmLocal_gCostFunction->get_Shifted_HGVolatilityParam()->get_ABCD();
//	general_info<<"a,b,c,d"<<std::endl;
//	general_info<<abcdParam.a_<<","<<abcdParam.b_<<","<<abcdParam.c_<<","<<abcdParam.d_<<","<<std::endl;
//	
//	// ------ print BAD SWAPTIONS
//	if(	!bad_swaptions_indices_.empty() )
//	{
//		general_info<<"BAD Calibrated ,,";
//		for(size_t i=0;i<bad_swaptions_indices_.size();++i)
//		{
//			general_info<<"("<<bad_swaptions_indices_[i].first<<";"<<bad_swaptions_indices_[i].second<<")  ,";
//		}
//		general_info<<std::endl;
//
//		general_info<<"BAD calibrated Value ,,";
//		for(size_t i=0;i<bad_swaptions_indices_.size();++i)
//		{
//			general_info<< calibrated_gDelegate_matrix_(bad_swaptions_indices_[i].first,bad_swaptions_indices_[i].second)<<" ,";
//		}
//	}
//	general_info<<std::endl;
//
//
//	if(isVirtualCalibration_)
//	{
//		general_info	<<max_param_rel_error_ <<", Max Param Error(%),"<<std::endl;
//	}
//
//	size_t nbRow = rows_stop_Error_.size();
//	assert(nbRow == rows_stop_Info_.size() );
//	general_info<<"Calib Row   ," ; for(size_t i=0;i<nbRow;++i) general_info<<(i+1)<<" , " ; general_info<<std::endl;
//	general_info<<"Stop Info   ," ; for(size_t i=0;i<nbRow;++i) general_info<<rows_stop_Info_[i]<<" , " ; general_info<<std::endl;
//	general_info<<"Calib Error ," ; for(size_t i=0;i<nbRow;++i) general_info<<rows_stop_Error_[i]<<" , " ; general_info<<std::endl;
//
//	base_general_result_info_ += general_info.str();
//}
//
//void CheyetteLocalCalibrator::printPlusPlus(const std::string& base_filename) const
//{
//	std::string path_OutPut = LMMPATH::get_output_path() + base_filename;
//
//	std::ofstream final_result ;
//
//	final_result.open(path_OutPut.c_str(), std::ios::app);
//
//	final_result<<"==============  General Results  ==============="<<std::endl;
//
//	final_result<<base_general_result_info_;
//
//	final_result<<std::endl;
//
//	final_result<<std::endl<<"==================== Quotation ERROR Matrix Results  ==================="<<std::endl;
//
//	final_result<<std::endl<<"## absolute error , "<<std::endl;
//	for(size_t i=1;i<abs_quote_error_matrix_.size1();++i)
//	{
//		for(size_t j=1;j<abs_quote_error_matrix_.size1()-i;++j)
//		{
//			final_result<<abs_quote_error_matrix_(i,j)<<",";
//		}
//		final_result<<std::endl;
//	}
//
//	final_result<<std::endl<<"## rel error (%), "<<std::endl;
//	for(size_t i=1;i<rel_quote_error_matrix_.size1();++i)
//	{
//		for(size_t j=1;j<rel_quote_error_matrix_.size1()-i;++j)
//		{
//			final_result<<rel_quote_error_matrix_(i,j)<<",";
//		}
//		final_result<<std::endl;
//	}
//



//void CheyetteDD_LocalCalibrator::minimizeNoConstraint(	QuantLib::Array xInitiate, 
//														QuantLib::LevenbergMarquardt minimizationSolver, 
//														CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR)
//{
//	QuantLib::Array					start_point = xInitiate ;
//	QuantLib::NoConstraint			noConstraint ;
//	QuantLib::Problem				optimizationProblem(*cheyetteBaseCostFunction_PTR, noConstraint, start_point); 
//	
//	QuantLib::EndCriteria::Type		endConvergenceType = minimizationSolver.minimize(optimizationProblem, stopCriteria_);	
//
//	std :: cout << " Criteria :" << endConvergenceType << std :: endl ;
//	std :: cout << " Root :" << optimizationProblem.currentValue() << std :: endl ;
//	std :: cout << " Min Function Value :" << optimizationProblem.functionValue () << std :: endl ;
//}

//calibratedArray : point de depart de l'algo d'optimisation
//calibratedArray est ensuite mis à jour avec la valeur calibrée
void CheyetteDD_LocalCalibrator::minimizePositiveConstraint(	QuantLib::Array& calibratedArray, 
																QuantLib::LevenbergMarquardt& minimizationSolver, 
																CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR)
{
	QuantLib::PositiveConstraint	positiveConstraint ;
	QuantLib::Problem				optimizationProblem(*cheyetteBaseCostFunction_PTR.get(), 
														positiveConstraint, 
														calibratedArray);  //pt de depart de l'algo
	
	QuantLib::EndCriteria::Type		endConvergenceType = minimizationSolver.minimize(optimizationProblem, stopCriteria_);	
	calibratedArray = optimizationProblem.currentValue() ;

	std::cout << " Criteria :" << endConvergenceType << std::endl ;
	std::cout << " Root :" << optimizationProblem.currentValue() << std::endl ;
	std::cout << " Min Function Value :" << optimizationProblem.functionValue () << std::endl ;
	std::cout << " --------------------------------------------------------------" << std::endl ;
	std::cout << std::endl ;
}



//calibration de sigma -> m -> sigma sur une swaption
void CheyetteDD_LocalCalibrator::calibrateOneSwaption(size_t indexSwaption)
{
	MarketData_PTR			marketData_PTR(cheyetteBaseCostFunction_->getMarketDataPTR()) ;
	CheyetteDD_VanillaSwaptionApproxPricer_PTR	cheyetteApprox_PTR(cheyetteBaseCostFunction_->getCheyetteDD_ApproxPricer_PTR()) ; 
	QuantLib::LevenbergMarquardt	minimizationSolver ; //(functionEpsilon(), rootEpsilon(), 1e-16);

	size_t nbMaxIterations = 2 ;
	for (size_t i = 0 ; i < nbMaxIterations ; ++i)
	{
		//minimisation sur sigma
		CheyetteDD_CostFunctionLevel_PTR	costFunctionSigma_PTR(new CheyetteDD_CostFunctionLevel(marketData_PTR, indexSwaption, cheyetteApprox_PTR)) ;
		minimizePositiveConstraint( calibrated_sigma_, minimizationSolver, costFunctionSigma_PTR) ;

		//minimisation sur m
		CheyetteDD_CostFunctionSkew_PTR		costFunctionSkew_PTR(new CheyetteDD_CostFunctionSkew(marketData_PTR, indexSwaption, cheyetteApprox_PTR)) ;	
		minimizePositiveConstraint( calibrated_m_, minimizationSolver, costFunctionSkew_PTR) ; 

		//minimisation sur sigma
		minimizePositiveConstraint( calibrated_sigma_, minimizationSolver, costFunctionSigma_PTR) ;
	}
	

//vol ATM modele
	double T			= cheyetteBaseCostFunction_->getMarketDataPTR()->get_aExpiry()[indexSwaption] ;
	double S0			= cheyetteApprox_PTR->get_buffer_s0_() ;
	double modelPrice	= cheyetteApprox_PTR->prixSwaptionApproxPiterbarg() ;
	double strike		= cheyetteBaseCostFunction_->getMarketDataPTR()->getStrikeATM()[indexSwaption] ;
	std::cout	<< "vol ATM : " << NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike, T) << std::endl ;

//vol ATM + shift modele 
	double shift = cheyetteBaseCostFunction_->getMarketDataPTR()->getShift() ;
	cheyetteApprox_PTR->setStrike(strike + shift) ;
	modelPrice	= cheyetteApprox_PTR->prixSwaptionApproxPiterbarg() ;
	std::cout	<< "vol ATM + shift : " << NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike+shift, T) << std::endl ;

//vol ATM + shift modele 
	cheyetteApprox_PTR->setStrike(strike - shift) ;
	modelPrice	= cheyetteApprox_PTR->prixSwaptionApproxPiterbarg() ;
	std::cout	<< "vol ATM - shift : " << NumericalMethods::Black_impliedVolatility(modelPrice, S0, strike-shift, T) << std::endl ;

	cheyetteApprox_PTR->setStrike(strike) ;
}

//boucle sur les differentes swaptions
void CheyetteDD_LocalCalibrator::calibrate()
{
	size_t nbSwaptions = cheyetteBaseCostFunction_->getMarketDataPTR()->getStrikeATM().size() ;
	std::vector<VanillaSwaption_PTR> vectSwaptions = cheyetteBaseCostFunction_->getMarketDataPTR()->getVect_swaptions() ;

	for (size_t indexSwaption = 0 ; indexSwaption < nbSwaptions ; ++indexSwaption)
	{
		//modifier l'approx pricer avec la swaption en cours
		cheyetteBaseCostFunction_->getCheyetteDD_ApproxPricer_PTR()->setSwaption(vectSwaptions[indexSwaption]) ;

		calibrateOneSwaption(indexSwaption) ; //calibration de sigma, puis de m et de sigma à nouveau
	}
	cheyetteBaseCostFunction_->getCheyetteDD_ApproxPricer_PTR()->get_CheyetteDD_Model()->show() ;



}

