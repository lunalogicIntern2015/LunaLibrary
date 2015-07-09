#include "TestCalibrator.h"

//void testCalib()
//{
////market Data
//	std::vector<size_t> a_expiry(1) ; a_expiry[0] = 1 ; //sur la tenor structure ?
//	std::vector<size_t> b_tenor(1) ; b_tenor[0] = 1 ;	
//	std::vector<double> strikeATM(1) ; strikeATM[0] =  2.78/100 ;  //1.2341/100 ; 0.36/100 ;
//	std::vector<double> volQuotes(1) ; volQuotes[0] =  29.41/100 ; //70.97/100 ;50.81/100 ;
//	std::vector<double> skew(1) ; skew[0] =  -7.4 ; //-3.8 ;
//	double shift = 5.0/10000 ;
//	CourbeInput_PTR courbe_PTR_test(createCourbeInput(5));
//
//	double strike = strikeATM[0] ;
//	LMM::Index  indexStart = 2 ; 
//	LMM::Index  indexEnd = 4 ; 
//	Tenor floatingLegTenorType = Tenor::_6M ;
//	Tenor fixedLegTenorType = Tenor::_1YR ;
//	LMMTenorStructure_PTR lmmTenorStructure(new LMMTenorStructure(Tenor::_6M, 30)) ;
//
//	VanillaSwap vanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, lmmTenorStructure); 
//
//	VanillaSwaption_PTR swaptionPTR(new VanillaSwaption(vanillaSwap , OptionType::OptionType::CALL));
//	std::vector<VanillaSwaption_PTR> vectSwaptions(1) ;
//	vectSwaptions[0] = swaptionPTR ;
//	MarketData_PTR mktData(new MarketData(a_expiry, b_tenor, strikeATM, volQuotes, skew, shift, vectSwaptions, "Black")) ;
//
//	
////Cheyette DD model
//
//	std::vector<double> x, sigma_y, m_y ;
//	x.push_back(0) ; x.push_back(1) ; 
//
//	double m_param = 0.4 ;
//	double k = 0.2 ;
//	double sigm = 0.20 ;
//
//	m_y.push_back(m_param) ; 
//	sigma_y.push_back(sigm) ;
//	
//	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
//	Piecewiseconst_RR_Function m(x, m_y) ; 
//
//	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
//	int shiftChoice = 1 ;
//	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
//	modele_test_PTR->show() ;
//
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR approx(
//					new CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, vectSwaptions[0]));
//
////calibrator
//	QuantLib::Size maxIterations	= 500 ;
//	QuantLib::Real rootEpsilon		= 1/10000 ;
//	QuantLib::Real functionEpsilon	= 1/10000 ;
//
//	size_t indexSwaption = 1 ;
//
//	CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR(
//				new CheyetteDD_CostFunctionLevel(mktData, indexSwaption, approx)) ;
//
//	Array sigmaInitiate(1) ; 
//	sigmaInitiate[0] = sigm ;
//	Array mInitiate(1) ;
//	mInitiate[0] = m_param ;
//	Array calibrated_sigma(1) ; 
//	calibrated_sigma = sigmaInitiate ;
//	Array calibrated_m(1) ;
//	calibrated_m = mInitiate ;
//
//	CheyetteDD_LocalCalibrator calibrator(	maxIterations,  rootEpsilon,   functionEpsilon,    
//											cheyetteBaseCostFunction_PTR,
//											sigmaInitiate, mInitiate,
//											calibrated_sigma, calibrated_m) ; 
//
//	calibrator.calibrate() ;
//
//}



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

//calibration
	marketData_Cheyette_LocalCalibration(config, pLmmSwaptionMarketData);

}


void marketData_Cheyette_LocalCalibration(	const CheyetteDD_CalibrationConfig& config,
											LmmSwaptionMarketData_PTR pLmmSwaptionMarketData)
{
	size_t nbYear = pLmmSwaptionMarketData->get_nbYear();
	std::string base_file_name = pLmmSwaptionMarketData->get_MarketDataBaseFileName();

	Tenor tenorfixedleg = Tenor::_1YR ;
	Tenor tenorfloatleg = Tenor::_6M  ;
	size_t fixedfloatRatio = tenorfixedleg.ratioTo(tenorfloatleg);

	std::string base_name;
	base_name = base_file_name+"_CheyetteDD_LocalCalibration" ;

	//create LMM components
	LMMTenorStructure_PTR pLMMTenorStructure( new LMMTenorStructure(tenorfloatleg,nbYear) );

//creation du Cheyette DD Model
	size_t shiftChoice = config.shiftChoice_ ;
	CourbeInput_PTR courbeInput = config.courbeInput_PTR_ ;

	CheyetteDD_Model_PTR cheyette_PTR(new CheyetteDD_Model(courbeInput, config.cheyetteDD_Param_, shiftChoice)) ;

//creation du ApproxPricer
	VanillaSwaption_PTR swaption_PTR(new VanillaSwaption()) ;
	CheyetteDD_VanillaSwaptionApproxPricer_PTR approxPricer(
									new CheyetteDD_VanillaSwaptionApproxPricer(cheyette_PTR, swaption_PTR));	
	
//creation CostFunction
	//pas un constructeur par defaut, transformer pLMM_market Data en CoTerminalQuotes

	CoTerminalSwaptionVol_CONSTPTR coTerminalSwaptionVol_PTR(new CoTerminalSwaptionVol()) ;
	CoTerminalSwaptionSkew_CONSTPTR coTerminalSwaptionSkew_PTR(new CoTerminalSwaptionSkew()) ;

	size_t indexSwaption = 0 ;
	CheyetteDD_CostFunctionLevel_PTR cheyetteCostFunctionLevel_PTR(
									new CheyetteDD_CostFunctionLevel(coTerminalSwaptionVol_PTR, indexSwaption, approxPricer) ) ;
	CheyetteDD_CostFunctionSkew_PTR cheyetteCostFunctionSkew_PTR(
									new CheyetteDD_CostFunctionSkew(coTerminalSwaptionSkew_PTR, indexSwaption, approxPricer) ) ;

// Create Calibrator + solve()
	size_t arraySize = config.cheyetteDD_Param_.sigma_.gety_().size() ;
	QuantLib::Array arraySigma(arraySize) ;
	QuantLib::Array arrayM(arraySize) ;
	CheyetteDD_LocalCalibrator cheyetteCalibrator(500, //maxIter
												 1e-8,   //x_epsilon (root epsilon)
												 1e-8,   //f_epsilon (function epsilon)   
												 arraySigma, arrayM, arraySigma, arrayM,
												 cheyetteCostFunctionLevel_PTR, cheyetteCostFunctionSkew_PTR) ;
	cheyetteCalibrator.solve() ;

//ecriture output
	std::ostringstream file_result_stream;file_result_stream<<base_name<<"_result.csv";
	std::string file_calibration_result(file_result_stream.str());
	//lmmCalibrator.printPlusPlus(file_calibration_result);
	
	//cheyetteCalibrator.printPlusPlus(file_calibration_result) ;

	//std::ostringstream file_gDelegate_stream;file_gDelegate_stream<<base_name<<"_gDelegate.csv";
	//std::string file_gDelegate_vol(file_gDelegate_stream.str() );
	//pGMatrixMapping->print(file_gDelegate_vol);

	//std::ostringstream file_vol_stream;file_vol_stream<<base_name<<"_vol.csv";
	//std::string file_calibrated_vol(file_vol_stream.str() );
	//pNoShifted_HGVolatilityParam->print( file_calibrated_vol );


	//{
	//	std::string common_result_file_name = "calib_result_CheyetteLocal.csv";

	//	std::string full_common_result_file = LMMPATH::get_Root_OutputPath() + common_result_file_name ;

	//	std::ofstream final_result ;
	//	final_result.open(full_common_result_file.c_str(), std::ios::app);

	//	final_result<<std::endl<<std::endl<< "============= Test At    "<<LMMPATH::get_TimeDateNow()
	//		<<",,,,,, Error LInf, "<<lmmCalibrator.get_QuoteError_LInf() <<std::endl ;
	//	final_result<< lmmCalibrator.get_BaseGeneral_Result_Info();

	//	final_result.close();	
	//}
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

CoTerminalSwaptionVol_PTR upperMatrixToCoTerminal(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData)
{

}