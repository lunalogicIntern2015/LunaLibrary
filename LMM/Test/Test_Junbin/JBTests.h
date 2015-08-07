#pragma once
#define BOOST_UBLAS_TYPE_CHECK 0
#include <ql/math/array.hpp>
#include <string>
#include <time.h>

#include <LMM/LmmSwaptionMarketData.h>
#include <LMM/calibration/LmmCalibrationConfig.h>
#include <LMM/Calibration/GMatrixMapping.h>
#include <LMM/Model/Shifted_HGVolatilityFunction.h>
#include <LMM/Pricer/McLmmPricer/McLmmVanillaSwapPricer.h>

#include <Instrument/Rate/Rate1.h>  
#include <Instrument/Rate/ConstRate.h>  
#include <Instrument/Rate/LiborRate.h>  
#include <Instrument/Rate/VanillaSwapRate.h>  
#include <LMM/Pricer/Longstaff_Schwartz/Basis.h>
#include <LMM/Pricer/Longstaff_Schwartz/Basis_Evaluator.h>
#include <LMM/Pricer/Longstaff_Schwartz/Regression_LS.h>
#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>
#include <Instrument/CallableOption/CallableInstrument.h>

#include <Instrument/InstrumentFactory.h>
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>
#include <LMM/Pricer/McLmmPricer/McLmmGenericSwapPricer.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_BackwardAlgo.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_ForwardAlgo.h>





//test JB
void vanillaSwapComparaisonExemple();
void Test_GeneticTargetSwapOneTrajectory();
void Test_McGeneticSwapLMMPricer();
void Test_McGeneticTargetSwapLMMPricing();
void JB_test_LmmCalibrationMarketData();
void Test_basis();
void Test_RegressionLS();
void Test_BermudanOption();
void Test_EV();
void Test_Longstaff_Schwartz_CallableSwap();
//void Test_LS();
void Test_Longstaff_Schwartz();
void Test_VanillaSwaptionPricing();
McLmm_PTR getMcLmmExample(LMMTenorStructure_PTR lmmTenorStructure, const std::vector<double>& initLiborValues, const LmmCalibrationConfig& config);
void Test_regression_LS();
void Test_pricing_forward_LS();
void Test_Longstaff_Schwartz_CallableSwap_for_test();
void Test_LS_pricing_allSubSet_basis();
Basis_CONSTPTR getSingleEVBasis_swaprateCall(const VanillaSwap& vanillaSwap, double strike, LMM::Index liborIndex);
Basis_CONSTPTR getSingleEVBasis_swaprate(const VanillaSwap& vanillaSwap, size_t swaprate_degree);
Basis_CONSTPTR getSingleEVBasis_Libor(LMM::Index liborIndex,  size_t liborRate_degree);
Basis_Evaluator_CONSTPTR getBasisEvaluator(const std::vector<size_t>& basis_evaluator_vect ,const McLmmVanillaSwapPricer& mcLmmVanillaSwapPricer);
Basis_CONSTPTR getBasis(const std::vector<size_t>& power_vect ,double coeff, const VanillaSwap& vanillaSwap, double strike, LMM::Index liborIndex);
void getAllSubsets(const vector<std::vector<size_t>>& set, std::vector<std::vector<std::vector<size_t>>>& subset);
void Test_LS_pricing_One_SubSet_basis(	const std::vector<std::vector<size_t>>& basis_subset, 
										double strike,
										LMM::Index  indexStart, 
										LMM::Index  indexEnd,
										const Tenor&	floatingLegTenorType,
										const Tenor&	fixedLegTenorType,
										const std::vector<double>& initLiborValues,
										const std::vector<LMM::Index>& exerciseDates,
										CallableInstrument_CONSTPTR callableGenericSwap,
										const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_backward,
										const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_forward,
										const std::vector<LMM::Index>& nbSimulation_vect,
										std::vector<std::vector<double>>& basis_value_on_allPath_buffer,
										ofstream& out,
										ofstream& out_test_time);
void representation_basis(const std::vector<size_t>& basis_vect, ofstream& out);
void representation_basis_subset(const std::vector<std::vector<size_t>>& basis_subset, ofstream& out);
void Test_10_basis();
void get_bisis_subset(std::vector<std::vector<std::vector<size_t>>>& subset);
void Test_LS_pricing_vol_correl(const LmmCalibrationConfig& config);
void Test_LS_pricing_parameter();

void Test_evaluation_basis();


typedef std::pair<Shifted_HGVolatilityFunction_PTR, GMatrixMapping_PTR>  GMatrix_Vol_gMapping;

LmmSwaptionMarketData_PTR JB_get_LmmSwaptionMarketData(const LmmCalibrationConfig& config, const std::string& input_file);
LMMTenorStructure_PTR JB_create___LMMTenorStructure_PTR(const size_t nbyear);
Correlation_PTR JB_create_InitCorrelation(const LmmCalibrationConfig& config);
Shifted_HGVolatilityFunction_PTR JB_marketData_LMM_ABCD_calibration(const LmmCalibrationConfig& config, LmmSwaptionMarketData_PTR pLmmSwaptionMarketData);
Correlation_PTR JB_marketData_LMM_Correlation_calibration(const LmmCalibrationConfig& config, LmmSwaptionMarketData_PTR pLmmSwaptionMarketData , const QuantLib::Array& found_abcd);
GMatrixMapping_PTR JB_marketData_LMM_CascadeExact_calibration( const LmmCalibrationConfig& config
														   , LmmSwaptionMarketData_PTR pLmmSwaptionMarketData 
														   , const QuantLib::Array& abcd_param 
														   , Correlation_PTR found_correlation_ptr 
														   );
GMatrix_Vol_gMapping JB_marketData_LMM_Global_gCalibration( const LmmCalibrationConfig& config
										, LmmSwaptionMarketData_PTR pLmmSwaptionMarketData 
										, Shifted_HGVolatilityFunction_PTR abcd_param 
										, Correlation_PTR found_correlation_ptr  
										, GMatrixMapping_PTR init_gMapping 
										);
void JB_marketData_LMM_Local_gCalibration( const LmmCalibrationConfig& config
									   , LmmSwaptionMarketData_PTR pLmmSwaptionMarketData 
									   , const QuantLib::Array& abcd_param 
									   , Correlation_PTR found_correlation_ptr  
									   , GMatrixMapping_PTR init_gMapping 
									   );

Shifted_HGVolatilityFunction_PTR		JB_marketData_LMM_shift_Calibration(	const LmmCalibrationConfig& config
																				, LmmSwaptionMarketData_PTR pLmmSwaptionMarketData 
																				, Shifted_HGVolatilityFunction_PTR param_h_g
																			);



//QuantLib::Array JB_marketData_LMM_ABCD_calibration_QL( const LmmCalibrationConfig& config
//													   , LmmSwaptionMarketData_PTR pLmmSwaptionMarketData
//													   );
//
//const Shifted_HGVolatilityParam::LowerTriangularMatrix& JB_marketData_LMM_Global_gCalibration( const LmmCalibrationConfig& config
//																	, LmmSwaptionMarketData_PTR pLmmSwaptionMarketData 
//																	, const QuantLib::Array& abcd_param
//																	, const Shifted_HGVolatilityParam::LowerTriangularMatrix& gMatrix
//																	, const std::vector<double>& shiftValues
//																	, Correlation_PTR found_correlation_ptr  
//																	, GMatrixMapping_PTR init_gMapping
//																	);
//
//
//const std::vector<double>&		JB_marketData_LMM_shift_Calibration_All(	const LmmCalibrationConfig& config
//																		, const QuantLib::Array& abcd_param
//																		, const Shifted_HGVolatilityParam::LowerTriangularMatrix& gMatrix
//																		, const std::vector<double>& shiftValues
//																		, LmmSwaptionMarketData_PTR pLmmSwaptionMarketData 
//																		);