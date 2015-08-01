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
void CheyetteDD_LocalCalibrator::minimizePositiveConstraint(QuantLib::Array& calibratedArray1D, 
															QuantLib::LevenbergMarquardt& minimizationSolver, 
															CheyetteDD_CostFunctionLevel pCostFunctionLevel)
{
	QuantLib::PositiveConstraint	positiveConstraint ;
	QuantLib::Problem				optimizationProblem(pCostFunctionLevel, 
														positiveConstraint, 
														calibratedArray1D);  //pt de depart de l'algo
	
	QuantLib::EndCriteria::Type		endConvergenceType = minimizationSolver.minimize(optimizationProblem, stopCriteria_);	
	calibratedArray1D = optimizationProblem.currentValue() ;

	std::cout << " Criteria :" << endConvergenceType << std::endl ;
	std::cout << " Root :" << optimizationProblem.currentValue() << std::endl ;
	std::cout << " Min Function Value :" << optimizationProblem.functionValue () << std::endl ;
	std::cout << " --------------------------------------------------------------" << std::endl ;

	o_ << " Criteria : ;" << endConvergenceType << std::endl ;
	o_ << " Root : ;" << optimizationProblem.currentValue() << std::endl ;
	o_ << " Min Function Value : ;" << optimizationProblem.functionValue () << std::endl ;
	o_ << std::endl ;
}

void CheyetteDD_LocalCalibrator::minimizeBoundaryConstraint(QuantLib::Array& calibratedArray1D, 
															QuantLib::LevenbergMarquardt& minimizationSolver, 
															CheyetteDD_CostFunctionSkew pCostFunctionSkew)
{
	//QuantLib::PositiveConstraint	bc ; //positiveConstraint ;
	QuantLib::BoundaryConstraint bc(0.000000001, 1.) ;
	//QuantLib::BoundaryConstraint bc(- 1., 1.) ;

	QuantLib::Problem				optimizationProblem(pCostFunctionSkew, 
														bc, 
														calibratedArray1D);  //pt de depart de l'algo
	
	QuantLib::EndCriteria::Type		endConvergenceType = minimizationSolver.minimize(optimizationProblem, stopCriteria_);	
	calibratedArray1D = optimizationProblem.currentValue() ;

	std::cout << " Criteria :" << endConvergenceType << std::endl ;
	std::cout << " Root :" << optimizationProblem.currentValue() << std::endl ;
	std::cout << " Min Function Value :" << optimizationProblem.functionValue () << std::endl ;
	std::cout << " --------------------------------------------------------------" << std::endl ;
	std::cout << std::endl ;

	o_ << " Criteria : ;" << endConvergenceType << std::endl ;
	o_ << " Root : ;" << optimizationProblem.currentValue() << std::endl ;
	o_ << " Min Function Value : ;" << optimizationProblem.functionValue () << std::endl ;
	o_ << std::endl ;
}

//solve le CheyetteDD_calibrator : minimisation 1D, une swaption
void CheyetteDD_LocalCalibrator::solve()
{

	//minimisations 1D sur sigma puis m
	QuantLib::LevenbergMarquardt	minimizationSolver ; //(functionEpsilon(), rootEpsilon(), 1e-16);

	//minimisation sur sigma
	minimizePositiveConstraint(calibrated_sigma_1D_, minimizationSolver, *costFunctionLevel_PTR_) ;  //.get()

	size_t nbMaxIterations = 2 ;
	for (size_t i = 0 ; i < nbMaxIterations ; ++i)
	{
		//minimisation sur m
		minimizeBoundaryConstraint(calibrated_m_1D_, minimizationSolver, *costFunctionSkew_PTR_) ; 

		//minimisation sur sigma
		minimizePositiveConstraint(calibrated_sigma_1D_, minimizationSolver, *costFunctionLevel_PTR_) ;
	}
	//mise à jour de calibrated_m_sigma et calibrated_m ok car passage par référence dans minimizeConstraint

//pour le smile
	double strike		= costFunctionLevel_PTR_->getCoTerminalSwaptionVol_PTR()->getStrike() ;
	double shift		= costFunctionSkew_PTR_->getCoTerminalSwaptionSkew_PTR()->getShift() ;

	o_ << " ; vol ATM - 5bp ; vol ATM ; vol ATM + 5bp " << std::endl ;
	o_ << "STRIKE ; " << strike - shift << " ; " << strike << " ; " << strike + shift << std::endl ;
//vol ATM - shift modele 
	double T			= costFunctionLevel_PTR_->getCoTerminalSwaptionVol_PTR()->getVectorExpiry() ;
	double S0			= costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->get_buffer_s0_() ;
	double annuity0		= costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->swapRateDenominator(0., 0., 0.) ;

	costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->setStrike(strike - shift) ;
	double modelPrice	= costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->prixSwaptionApproxPiterbarg() ;
	double volATMmoins	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
																			S0, strike - shift, T) ;
	o_ << "VOL MODEL ; " << volATMmoins ;

//vol ATM 
	costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->setStrike(strike) ;
	modelPrice	= costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->prixSwaptionApproxPiterbarg() ;
	double volATM	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
																			S0, strike, T) ;
	o_ << " ; " << volATM ;

//vol ATM + shift modele 
	costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->setStrike(strike + shift) ;
	modelPrice	= costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->prixSwaptionApproxPiterbarg() ;
	double volATMplus	= NumericalMethods::Black_SwaptionImpliedVolatility(modelPrice, annuity0,   							  
																			S0, strike + shift, T) ;
	o_ << " ; " << volATMplus << std::endl ;


	costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->setStrike(strike) ;

	//costFunctionLevel_PTR_->getCheyetteDD_ApproxPricer_PTR()->get_CheyetteDD_Model()->show() ;
	//retrieve_calib_info();
}

//void CheyetteDD_LocalCalibrator::printPlusPlus(const std::string& base_filename) const
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
//
//	///////  PRINT model quote and market quotes
//	UpperTriangularDoubleMatrix mkt_black_vol_matrix = pLmmBaseCostFunction_->get_MarketQuoteValuesMatrix();
//	UpperTriangularDoubleMatrix mdl_black_vol_matrix = pLmmBaseCostFunction_->get_ModelQuoteValuesMatrix();
//
//	final_result<<std::endl<<"## Model Black Vol, "<<std::endl;
//	for(size_t i=1;i<mdl_black_vol_matrix.size1();++i)
//	{
//		for(size_t j=1;j<mdl_black_vol_matrix.size1()-i;++j)
//		{
//			final_result<<mdl_black_vol_matrix(i,j)<<",";
//		}
//		final_result<<std::endl;
//	}
//
//	final_result<<std::endl<<"## Market Black Vol, "<<std::endl;
//	for(size_t i=1;i<mkt_black_vol_matrix.size1();++i)
//	{
//		for(size_t j=1;j<mkt_black_vol_matrix.size1()-i;++j)
//		{
//			final_result<<mkt_black_vol_matrix(i,j)<<",";
//		}
//		final_result<<std::endl;
//	}
//
//	final_result<<std::endl<<"==================== Quotation ERROR Vector Results  ==================="<<std::endl;
//	final_result<<"(Expir;Tenor) , , MKT Quotes, MDL Quotes, , Abs Error, Rel Error (%)"<<std::endl;
//
//
//	for(size_t iExpirity = 1; iExpirity<mkt_black_vol_matrix.size1()-1; ++iExpirity) // row
//	{
//		size_t nb_jSwaption = mkt_black_vol_matrix.size2()- iExpirity;
//		for(size_t jTenor = 1;jTenor< nb_jSwaption; ++jTenor)
//		{
//			final_result<<"("<<iExpirity<<";"<<jTenor<<"), ,";
//			final_result<<mkt_black_vol_matrix(iExpirity,jTenor) <<","
//				<<mdl_black_vol_matrix(iExpirity,jTenor) <<", ,"
//				<<abs_quote_error_matrix_(iExpirity,jTenor) <<","
//				<<rel_quote_error_matrix_(iExpirity,jTenor) <<",";
//
//			if(std::abs(rel_quote_error_matrix_(iExpirity,jTenor))>30)
//			{
//				final_result<<", BAD BAD,"<<std::endl;
//			}
//			else if(std::abs(rel_quote_error_matrix_(iExpirity,jTenor))>10)
//			{
//				final_result<<", BAD,"<<std::endl;
//			}
//			else
//			{
//				final_result<<","<<std::endl;
//			}			
//		}
//	}
//	final_result<<std::endl
//		<<", , , , Max Error ,"<< max_quote_abs_error_ <<","<< max_quote_rel_error_        <<","<<std::endl<<std::endl;
//
//	if(isVirtualCalibration_)
//	{
//		final_result<<std::endl<<"==================== Param ERROR Vector Results for Virtual Test gDelegate ==================="<<std::endl;
//		final_result<<"(Expir;Tenor) , , TrueParam, Calibrated Param, Error (%)"<<std::endl;
//
//
//		for(size_t iExpirity = 1; iExpirity<mkt_black_vol_matrix.size1()-1; ++iExpirity) // row
//		{
//			size_t nb_jSwaption = mkt_black_vol_matrix.size2()- iExpirity;
//			for(size_t jTenor = 1;jTenor< nb_jSwaption; ++jTenor)
//			{
//				final_result<<"("<<iExpirity<<";"<<jTenor<<"), ,";
//				final_result<<true_gDelegate_matrix_(iExpirity,jTenor) <<","
//					<<calibrated_gDelegate_matrix_(iExpirity,jTenor) <<", ,"
//					<<rel_error_gDelegate_matrix_(iExpirity,jTenor) <<",";
//
//				if(std::abs(rel_error_gDelegate_matrix_(iExpirity,jTenor))>15)
//				{
//					final_result<<", BAD BAD,"<<std::endl;
//				}
//				else if(std::abs(rel_error_gDelegate_matrix_(iExpirity,jTenor))>5)
//				{
//					final_result<<", BAD,"<<std::endl;
//				}
//				else
//				{
//					final_result<<","<<std::endl;
//				}			
//			}
//		}
//		final_result<<std::endl
//			<<", , , ,Max Error(%) ,"<< max_param_rel_error_        <<","<<std::endl<<std::endl;
//	}
//
//	final_result.close();
//}
