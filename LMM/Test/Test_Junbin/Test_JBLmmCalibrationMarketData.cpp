#include <LMM/Test/Test_Junbin/JBTests.h>
#include <LMM/Test/Test_CalibrationConfig.h>

#include <iostream>
#include <cassert>
#include <string.h>
#include <cmath>
#include <fstream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <ctime>
#include <cmath>

// ---- include for QuantLib calibration -------
#include <ql/termstructures/volatility/abcdcalibration.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/bfgs.hpp> 
#include <ql/math/optimization/conjugategradient.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
// ---- include for QuantLib calibration -------

#include <LMM/Helper/GenericPath.h>
#include <LMM/Helper/TenorType.h>
#include <LMM/Helper/LMMTenorStructure.h>
#include <Numeric/NumericalMethods.h>
#include <RNGenerator/McGenerator.h>
#include <LMM/Helper/Noise.h>
#include <LMM/Helper/InputFileManager.h>
#include <LMM/Helper/LmmGnuplotPrinterMatrix.h>

#include <LMM/calibration/ATMSwaptionMarketData.h>
#include <LMM/calibration/SwaptionMarketDataContainer.h>
#include <LMM/calibration/SwaptionMarketDataManager.h>

#include <LMM/Model/Lmm.h>
#include <LMM/Model/Correlation.h>
#include <LMM/Model/ConstShifted_HGVolatilityFunction.h>
#include <LMM/Calibration/Calibrator_CostFunction/LmmBaseCalibrator.h>
#include <LMM/Calibration/Calibrator_CostFunction/LmmABCDCalibrator.h>
#include <LMM/Calibration/Calibrator_CostFunction/LmmSkewCostFunction.h> 
#include <LMM/Pricer/LmmApproximationPricer/LmmVanillaSwaptionApproxPricer_Rebonato.h>

#include <LMM/LmmSwaptionMarketData.h>

#include <LMM/Calibration/Calibrator_CostFunction/LmmShiftCalibrator.h>

#include <LMM/Calibration/Calibrator_CostFunction/LmmGlobal_gCalibrator.h>
#include <LMM/Calibration/Calibrator_CostFunction/LmmGlobal_gCostFunction.h>

#include <LMM/Model/Shifted_HGVolatilityFunction.h>





void JB_test_calib_gMatrix_OneFile(const std::string& mkt_data_file);


void JB_test_calib_shift_gMatrix_OneFile(const std::string& mkt_data_file);

void JB_test_LmmCalibrationMarketData()
{

	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
	JB_test_calib_shift_gMatrix_OneFile(mkt_file_list[0]);

	//////series of tests all calibator, with all data files, and allowing negative values of g .
	//JB_test_calib_gMatrixNegative_allData();
}



void JB_test_calib_shift_gMatrix_OneFile(const std::string& mkt_data_file)
{

	//get input file 
	std::string folder_name;// = "TotalCalib\\" ;  config.use_positive_constraint_=true;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );



	LmmCalibrationConfig config;  //JB: default parameters set: model_nbYear_=16, correl_FullRank_ = 2*model_nbYear_ + 1  
									// case fixed/float = 2 typically lmm 6MO/1YR

	config.model_nbYear_=16;		// set model_nbYear_
	config.correl_FullRank_=2*config.model_nbYear_+1;

	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData=JB_get_LmmSwaptionMarketData(config, mkt_data_file);

	config.correl_ReducedRank_= 3 ;
	config.correl_alpha_ = 0.000000001 ; 
	config.correl_beta_  = 0.05;
	Shifted_HGVolatilityFunction_PTR  shifted_HGVolatilityFunction= JB_marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);

	//Correlation_PTR found_correlation_ptr = marketData_LMM_Correlation_calibration(config,pLmmSwaptionMarketData,found_abcd);

	

	//marketData_LMM_CascadeExact_calibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config) );
	//config.use_local_calib_=true;
	//marketData_LMM_Local_gCalibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config), GMatrixMapping_PTR() );
	//config.use_local_calib_=false;
	config.penalty_time_homogeneity_ = 1e-4 ; 
	config.penalty_libor_ = 1e-6 ; 
	config.use_positive_constraint_= true;

	GMatrix_Vol_gMapping gMatrix_Vol_gMapping	=	JB_marketData_LMM_Global_gCalibration(config, pLmmSwaptionMarketData, shifted_HGVolatilityFunction, create_InitCorrelation(config), GMatrixMapping_PTR() );
	
	Shifted_HGVolatilityFunction_PTR param_h_g_shifted	=	JB_marketData_LMM_shift_Calibration(	config
																								, pLmmSwaptionMarketData 
																								, gMatrix_Vol_gMapping.first
																								);

	JB_marketData_LMM_Global_gCalibration(config, pLmmSwaptionMarketData, param_h_g_shifted , create_InitCorrelation(config), gMatrix_Vol_gMapping.second);
	cout << "end of program" << endl;
//		JB_marketData_LMM_Global_gCalibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config), GMatrixMapping_PTR() )
}

