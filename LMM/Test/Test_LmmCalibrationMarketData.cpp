#include <LMM/Test/Tests.h>
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

#include <LMM/helper/GenericPath.h>
#include <LMM/helper/TenorType.h>
#include <LMM/helper/LMMTenorStructure.h>
#include <LMM/numeric/NumericalMethods.h>
#include <LMM/RNGenerator/McGenerator.h>
#include <LMM/helper/Noise.h>
#include <LMM/helper/InputFileManager.h>
#include <LMM/helper/LmmGnuplotPrinterMatrix.h>

#include <LMM/calibration/ATMSwaptionMarketData.h>
#include <LMM/calibration/SwaptionMarketDataContainer.h>
#include <LMM/calibration/SwaptionMarketDataManager.h>

#include <LMM/LmmModel/Lmm.h>
#include <LMM/LmmModel/Correlation.h>
#include <LMM/LmmModel/ConstShifted_HGVolatilityFunction.h>
#include <LMM/LmmModel/LmmABCDCostFunction.h>
#include <LMM/LmmModel/LmmABCDCalibrator.h>
#include <LMM/LmmModel/LmmSkewCostFunction.h>
#include <LMM/pricer/LmmVanillaSwaptionApproxPricer_Rebonato.h>

#include <LMM/LmmModel/LmmSwaptionMarketData.h>

#include <LMM/LmmModel/LmmShiftCalibrator.h>

#include <LMM/LmmModel/LmmGlobal_gCalibrator.h>
#include <LMM/LmmModel/LmmGlobal_gCostFunction.h>
#include <LMM/LmmModel/LmmLocal_gCalibrator.h>
#include <LMM/LmmModel/LmmLocal_gCostFunction.h>
#include <LMM/LmmModel/LmmCascade_gCalibrator.h>
#include <LMM/LmmModel/LmmCascade_gCostFunction.h>

//pour Cheyette
#include <Cheyette/Calibrator/CheyetteDD_CalibrationConfig.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Pricer/CheyetteDD_VanillaSwaptionApproxPricer.h>
#include <Cheyette/Calibrator/CheyetteBaseCalibrator.h>

void test_calib_gMatrix_OneFile(const std::string& mkt_data_file);

void test_calib_gMatrixNegative_allData();

void test_LmmCalibrationMarketData()
{

	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
	test_calib_gMatrix_OneFile(mkt_file_list[0]);

	//////series of tests all calibator, with all data files, and allowing negative values of g .
	//test_calib_gMatrixNegative_allData();
}


void test_calib_gMatrix_OneFile(const std::string& mkt_data_file)
{
	std::string folder_name;// = "TotalCalib\\" ;  config.use_positive_constraint_=true;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	LmmCalibrationConfig config;
	config.model_nbYear_=16;

	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData=get_LmmSwaptionMarketData(config, mkt_data_file);

	config.correl_ReducedRank_= 30 ; config.correl_alpha_ = 0.000000001 ; config.correl_beta_  = 0.05;
	QuantLib::Array found_abcd = marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);

	//Correlation_PTR found_correlation_ptr = marketData_LMM_Correlation_calibration(config,pLmmSwaptionMarketData,found_abcd);

	config.correl_ReducedRank_= 3 ;

	//marketData_LMM_CascadeExact_calibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config) );
	//config.use_local_calib_=true;
	//marketData_LMM_Local_gCalibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config), GMatrixMapping_PTR() );
	//config.use_local_calib_=false;
	config.penalty_time_homogeneity_ = 1e-4 ; config.penalty_libor_ = 1e-6 ; config.use_positive_constraint_= true;
	marketData_LMM_Global_gCalibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config), GMatrixMapping_PTR() );
}


void test_calib_gMatrixNegative_allData()
{
	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
	//size_t nbFile = 2;
	size_t nbFile = mkt_file_list.size();

	LmmCalibrationConfig config;
	config.model_nbYear_=16;

	//std::cout<<"------------------------------------------------------  HORIZON "<<nbYear<<" YR "<<std::endl;//ctndebug_calib_todelete
	for(size_t iFile=0 ; iFile< nbFile ; ++iFile)
	{
		std::string folder_name;// = "TotalCalib\\" ;  config.use_positive_constraint_=true;
		std::string base_name_file = LMMPATH::get_BaseFileName(mkt_file_list[iFile]) + "\\";
		folder_name+=base_name_file;
		LMMPATH::reset_Output_SubFolder(folder_name );

		LmmSwaptionMarketData_PTR pLmmSwaptionMarketData=get_LmmSwaptionMarketData(config, mkt_file_list[iFile]);

		config.correl_ReducedRank_= 30 ; config.correl_alpha_ = 0.000000001 ; config.correl_beta_  = 0.05;
		QuantLib::Array found_abcd = marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);

		//Correlation_PTR found_correlation_ptr = marketData_LMM_Correlation_calibration(config,pLmmSwaptionMarketData,found_abcd);

		config.correl_ReducedRank_= 3 ;

		//marketData_LMM_CascadeExact_calibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config) );
		//config.use_local_calib_=true;
		//marketData_LMM_Local_gCalibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config), GMatrixMapping_PTR() );
		//config.use_local_calib_=false;
		//config.penalty_time_homogeneity_ = 1e-4 ; config.penalty_libor_ = 1e-6 ; config.use_positive_constraint_=false;
		marketData_LMM_Global_gCalibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config), GMatrixMapping_PTR() );

	}
}


/* Cheyette */

void test_CheyetteCalibrationMarketData()
{
	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
	test_calib_Cheyette_OneFile(mkt_file_list[0]);
}


void test_calib_Cheyette_OneFile(const std::string& mkt_data_file)
{
	std::string folder_name;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	CheyetteDD_CalibrationConfig config;
	config.model_nbYear_=16;

//recuperation et stockage des donnees de marche OK pour Cheyette
	//convertir la UpperTrianlge Matrix en Coterminal SwaptionQuote
	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketData(config.model_nbYear_, mkt_data_file);

//transformation sous Array des parametres de CheyetteDD
	QuantLib::Array found_abcd = marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);

	marketData_Cheyette_LocalCalibration(config, pLmmSwaptionMarketData, found_abcd);

//	config.use_local_calib_=true;
//	config.use_positive_constraint_= true;

}


void marketData_Cheyette_LocalCalibration(	const CheyetteDD_CalibrationConfig& config,
											LmmSwaptionMarketData_PTR pLmmSwaptionMarketData,
											const QuantLib::Array& arraySigma, 
											const QuantLib::Array& arrayM)
{
	//assert(config.use_local_calib_);
	size_t nbYear = pLmmSwaptionMarketData->get_nbYear();
	std::string base_file_name = pLmmSwaptionMarketData->get_MarketDataBaseFileName();

	Tenor tenorfixedleg = Tenor::_1YR ;
	Tenor tenorfloatleg = Tenor::_6M  ;
	size_t fixedfloatRatio = tenorfixedleg.ratioTo(tenorfloatleg);

	std::string base_name;
	base_name = base_file_name+"_CheyetteDD_LocalCalibration" ;

	//create LMM components
	LMMTenorStructure_PTR pLMMTenorStructure( new LMMTenorStructure(tenorfloatleg,nbYear) );


//creation du struct param : 
	//const double a=abcd_param[0];
	//const double b=abcd_param[1];
	//const double c=abcd_param[2];
	//const double d=abcd_param[3];
	//Shifted_HGVolatilityParam::ABCDParameter abcdParam(a,b,c,d);
	//ConstShifted_HGVolatilityParam_PTR pNoShifted_HGVolatilityParam( 
	//	new ConstShifted_HGVolatilityParam(pLMMTenorStructure, abcdParam, 1., 0.));

		/*	CheyetteDD_Parameter(double k, 
			                 const Piecewiseconst_RR_Function& sigma, 
							 const Piecewiseconst_RR_Function& m)		*/

	std::vector<double> vectSigma	= arrayToVector(arraySigma) ;
	std::vector<double> vectM		= arrayToVector(arrayM) ;
	std::vector<double> x			= config.cheyetteDD_Param_.sigma_.getx_() ;
	double k = config.cheyetteDD_Param_.k_ ;
	Piecewiseconst_RR_Function sigmaFunc(x, vectSigma) ;
	Piecewiseconst_RR_Function mFunc(x, vectM) ;
	CheyetteDD_Model::CheyetteDD_Parameter cheyetteParam(k, sigmaFunc, mFunc) ;

//creation du modele LMMptr / Cheyette DD Model
//	Lmm_PTR lmm_ptr(new Lmm(dispersion) );
	size_t shiftChoice = config.shiftChoice_ ;
	CourbeInput_PTR courbeInput = config.courbeInput_PTR_ ;

	CheyetteDD_Model_PTR cheyette_PTR(new CheyetteDD_Model(courbeInput, cheyetteParam, shiftChoice)) ;

//creation du ApproxPricer
	CheyetteDD_VanillaSwaptionApproxPricer_PTR approxPricer(
									new CheyetteDD_VanillaSwaptionApproxPricer(cheyette_PTR, ));	
	
//creation CostFunction
	//LmmBaseCostFunction_PTR pLmmCostFunction(new LmmLocal_gCostFunction
	//	(
	//	pLmmVanillaSwaptionApproxPricer_Rebonato, 
	//	pLmmSwaptionMarketData->get_LiborQuotes(), 
		pLmmSwaptionMarketData->get_SwaptionQuotes_ATM(),
	//	) );

	CheyetteBaseCostFunction_PTR cheyetteCostFunction_PTR(
									new CheyetteBaseCostFunction( , , approxPricer) ) ;
		//CheyetteBaseCostFunction(	MarketData_PTR marketData_PTR, size_t indexSwaption,  
		//							CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR)

//costumize swaptions weights
	////UpperTriangularDoubleMatrix swpm_weight_matrix = pLmmCostFunction->get_SwaptionWeightMatrix();
	////swpm_weight_matrix(7,1)=1e-6;
	////swpm_weight_matrix(10,1)=1e-6;
	////swpm_weight_matrix(5,3)=0.;
	//pLmmCostFunction->reset_SwaptionWeightMatrix(swpm_weight_matrix);

// Create Calibrator
	LmmLocal_gCalibrator lmmCalibrator
		(
		*pGMatrixMapping.get()
		, 3000 //maxIter
		, 1e-11   //x_epsilon
		, 1e-11   //f_epsilon    
		, pLmmCostFunction
		);

	if(config.use_positive_constraint_)
		lmmCalibrator.activate_PositiveConstraint();

	lmmCalibrator.solve();

	std::ostringstream file_result_stream;file_result_stream<<base_name<<"_result.csv";
	std::string file_calibration_result(file_result_stream.str());
	lmmCalibrator.printPlusPlus(file_calibration_result);

	std::ostringstream file_gDelegate_stream;file_gDelegate_stream<<base_name<<"_gDelegate.csv";
	std::string file_gDelegate_vol(file_gDelegate_stream.str() );
	pGMatrixMapping->print(file_gDelegate_vol);

	std::ostringstream file_vol_stream;file_vol_stream<<base_name<<"_vol.csv";
	std::string file_calibrated_vol(file_vol_stream.str() );
	pNoShifted_HGVolatilityParam->print( file_calibrated_vol );


	{
		std::string common_result_file_name = "calib_result_gLocal.csv";

		std::string full_common_result_file = LMMPATH::get_Root_OutputPath() + common_result_file_name ;

		std::ofstream final_result ;
		final_result.open(full_common_result_file.c_str(), std::ios::app);

		final_result<<std::endl<<std::endl<< "============= Test At    "<<LMMPATH::get_TimeDateNow()
			<<",,,,,, Error LInf, "<<lmmCalibrator.get_QuoteError_LInf() <<std::endl ;
		final_result<< lmmCalibrator.get_BaseGeneral_Result_Info();

		final_result.close();	
	}
}


std::vector<double> arrayToVector(QuantLib::Array a)
{
	std::vector<double> vect ;
	vect.reserve(a.size()) ;
	for (size_t i = 0 ; i < a.size() ; ++i)
	{
		vect[i] = a[i] ;
	}
	return vect ;
}

