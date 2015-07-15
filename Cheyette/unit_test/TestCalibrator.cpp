#include "TestCalibrator.h"
#include <cassert>

//data converter

std::vector<size_t> getVectorExpiry(size_t coterminal)
{
	std::vector<size_t>	vectExpiry(coterminal) ;
	vectExpiry[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		vectExpiry[i] = i ;
	}
	return vectExpiry ;
}
std::vector<size_t> getVectorTenor(size_t coterminal)
{
	std::vector<size_t>	vectTenor(coterminal) ;
	vectTenor[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		vectTenor[i] = coterminal - i ;
	}
	return vectTenor ;
}

std::vector<double> getVectorStrike(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal)
{
	std::vector<double>	vectStrike(coterminal) ;
	vectStrike[0] = 0 ;
	UpperTriangularDoubleMatrix strike = pLmmSwaptionMarketData->get_SwaptionQuotes_ATM()->get_UpperTriangularStrike() ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		vectStrike[i] = strike(i, coterminal - i) ;
	}
	return vectStrike ;
}


std::vector<double> getVectorVol(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal)
{
	// 1st row and column not used! size = nbYear + 1
	size_t nbLignes		= pLmmSwaptionMarketData->get_SwaptionQuotes_ATM()->size1() ;	//nb year +1
	assert(coterminal <= nbLignes) ;

	std::vector<double>	vectVol(coterminal) ;
	vectVol[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		//boost::numeric::ublas::matrix<SwaptionQuote >
		UpperTriangularVanillaSwaptionQuotes quotes = 
					pLmmSwaptionMarketData->get_SwaptionQuotes_ATM()->get_UpperTriangularVanillaSwaptionQuotes() ;
		SwaptionQuote qu = quotes(i, coterminal - i) ;			//typedef std::pair<VanillaSwaption, double> SwaptionQuote; 
		vectVol[i] = qu.second ; 
	}
	return vectVol ;
}

std::vector<double> getVectorSkew(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal)
{
	// 1st row and column not used! size = nbYear + 1
	size_t nbLignes		= pLmmSwaptionMarketData->get_SwaptionQuotes_skew()->size1() ;	//nb year +1
	assert(coterminal <= nbLignes) ;

	std::vector<double>	vectSkew(coterminal) ;
	vectSkew[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		//boost::numeric::ublas::matrix<SwaptionQuote >
		UpperTriangularVanillaSwaptionQuotes quotes = 
					pLmmSwaptionMarketData->get_SwaptionQuotes_skew()->get_UpperTriangularVanillaSwaptionQuotes() ;
		SwaptionQuote qu = quotes(i, coterminal - i) ;			//typedef std::pair<VanillaSwaption, double> SwaptionQuote; 
		vectSkew[i] = qu.second ; 
	}
	return vectSkew ;
}



CourbeInput upperMatrixZC_To_CourbeInput_yield(LiborQuotes_PTR liborQuotes_PTR)
{
	std::vector<double> listeMatu(liborQuotes_PTR->get_LMMTenorStructure_PTR()->get_tenorDate()) ; 	
	std::vector<double> ZC(liborQuotes_PTR->get_VectDiscountFactor()) ;
	std::vector<double> tauxZC(ZC.size()) ;

	for (size_t i = 1 ; i < tauxZC.size() ; ++i)
	{
		tauxZC[i] = - 1./listeMatu[i] * log(ZC[i]) ;
	}
	tauxZC[0] = tauxZC[1] ;
	return CourbeInput(listeMatu, tauxZC);
}

std::vector<VanillaSwaption_PTR> getVectorSwaptions(const Tenor& tenorFloat, const Tenor& tenorFixed,
												const std::vector<size_t>& vectExpiry, 
												const std::vector<size_t>& vectTenor, 
												const std::vector<double>& vectStrike)
{
	std::vector<VanillaSwaption_PTR> vectSwaptions(vectExpiry.size()) ;
	VanillaSwaption_PTR v(new VanillaSwaption()) ;
	vectSwaptions[0] = v ;
	for (size_t i = 1 ; i < vectSwaptions.size() ; ++i)
	{
		size_t indexStart					= vectExpiry[i] / tenorFloat.YearFraction() ;
		size_t indexEnd						= indexStart + vectTenor[i] / tenorFloat.YearFraction() ;
		LMMTenorStructure_PTR pTenorStructure(new LMMTenorStructure(tenorFloat, 40)) ;
		VanillaSwap swap = VanillaSwap(vectStrike[i], indexStart, indexEnd, tenorFloat, tenorFixed, pTenorStructure) ;
		VanillaSwaption_PTR v(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
		vectSwaptions[i] = v ;				
	}
	return vectSwaptions ;
}

//swaption aYbY, coterminal = a + b
//calibre sur toutes les swaptions coterminales
void testCalib2(size_t fileNumber, size_t coterminal)
{
//lecture du fichier et recuperation des market data
	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();

	std::string mkt_data_file = mkt_file_list[fileNumber] ;

	std::string folder_name;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	size_t model_nbYear = 16 ; //definit jusqu'où est lu VCUB
	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketData(model_nbYear, mkt_data_file);

//recuperation de la courbe taux, tenor, expiry, vol, skew, swaptions (vector) à partir de LmmSwaptionMarketData_PTR
	std::vector<size_t>	vectExpiry = getVectorExpiry(coterminal) ;
	std::vector<size_t>	vectTenor = getVectorTenor(coterminal) ;
	
	std::vector<double>	vectStrike = getVectorStrike(pLmmSwaptionMarketData, coterminal) ;
	std::vector<double>	vectVol = getVectorVol(pLmmSwaptionMarketData, coterminal) ;	
	std::vector<double>	vectSkew = getVectorSkew(pLmmSwaptionMarketData, coterminal) ;		
	
	CourbeInput_PTR courbeInput_PTR(new 
		CourbeInput(upperMatrixZC_To_CourbeInput_yield(pLmmSwaptionMarketData->get_LiborQuotes()) )
								) ;
	
	Tenor tenorFloat = Tenor::_6M ;
	Tenor tenorFixed = Tenor::_12M ;
	std::vector<VanillaSwaption_PTR> vectSwaptions = getVectorSwaptions(tenorFloat, tenorFixed,
																	vectExpiry, vectTenor, vectStrike) ;

	double shift = 5.0/10000 ;
	MarketData_PTR mktData(new MarketData(vectExpiry, vectTenor, vectStrike, vectVol, vectSkew, shift, vectSwaptions, "Black"));
	

//Cheyette DD model
	double k = 0.2 ;
	double sigma = 0.20 ;
	double m = 0.4 ;	

	Piecewiseconst_RR_Function sigmaFunc(coterminal - 1, sigma) ; 
	Piecewiseconst_RR_Function mFunc(coterminal - 1, m) ; 

	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigmaFunc, mFunc) ;
	int shiftChoice = 3 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbeInput_PTR, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//ecriture dans fichier
	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_output_path() ;
	fileName_s << directory << "calib.csv" ; 
	std::string fileName = fileName_s.str();

	ofstream o;

	o.open(fileName,  ios::out | ios::app );
	o	<<	endl;
	modele_test_PTR->print(o) ;
	switch (shiftChoice)
	{
		case 1 :
			o << "shift 1 : r(t) / r(0)" << endl ;
			break ;
		case 2 :
			o << "shift 2 : r(t) / f(0, t)" << endl ;
			break ;
		case 3 :
			o << "shift 3 : S(t) / S(0)" << endl ;
			break ;
		default:
			std::cout << "invalide shiftChoice" << std::endl ;
			throw "exception" ;
			break;
	}
	
	courbeInput_PTR->print(o) ;
	o << endl ;

//calibrator
	QuantLib::Size maxIterations	= 500 ;
	QuantLib::Real rootEpsilon		= 1./10000 ;
	QuantLib::Real functionEpsilon	= 1./10000 ;

	//minimisation 1D, on passe en parametre sigma[indexSwaption] ou m[indexSwaption]
	for (size_t indexSwaption = 1 ; indexSwaption < vectSwaptions.size() ; ++indexSwaption)
	{
		CheyetteDD_VanillaSwaptionApproxPricer_PTR approx(
					new CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, vectSwaptions[indexSwaption]));

		CoTerminalSwaptionVol_CONSTPTR coterminalSwaptionVol_PTR 
				(new CoTerminalSwaptionVol(	vectVol[indexSwaption], 
											vectExpiry[indexSwaption], 
											vectTenor[indexSwaption], 
											vectStrike[indexSwaption], 
											vectSwaptions[indexSwaption])) ;
		CheyetteDD_CostFunctionLevel_PTR cheyetteBaseCostFunctionLevel_PTR(
					new CheyetteDD_CostFunctionLevel(o, coterminalSwaptionVol_PTR, approx, indexSwaption)) ;

		CoTerminalSwaptionSkew_CONSTPTR coTerminalSwaptionSkew_PTR
				(new CoTerminalSwaptionSkew(vectSkew[indexSwaption], 
											vectExpiry[indexSwaption], 
											vectTenor[indexSwaption], 
											vectStrike[indexSwaption], 
											vectSwaptions[indexSwaption], 
											shift)) ;
		CheyetteDD_CostFunctionSkew_PTR cheyetteBaseCostFunctionSkew_PTR(
					new CheyetteDD_CostFunctionSkew(o, coTerminalSwaptionSkew_PTR, approx, indexSwaption)) ;

		std::vector<double> v_sigma(1) ; v_sigma[0] = sigma ;
		Array sigmaInitiate = vectorToArray(v_sigma) ; 

		std::vector<double> v_m(1) ; v_m[0] = m ;
		Array mInitiate = vectorToArray(v_m) ;

		Array calibrated_sigma(sigmaInitiate) ;
		Array calibrated_m(mInitiate) ;

		CheyetteDD_LocalCalibrator calibrator(	o, 
												maxIterations,  rootEpsilon,   functionEpsilon,    
												sigmaInitiate, mInitiate,
												calibrated_sigma, calibrated_m,
												cheyetteBaseCostFunctionLevel_PTR,
												cheyetteBaseCostFunctionSkew_PTR) ; 

		calibrator.solve() ;	
	}
	modele_test_PTR->print(o) ;
	o.close() ;
}


//void testCalib()
//{
////market Data
//	std::vector<size_t> a_expiry(1) ; a_expiry[0] = 1 ; //sur la tenor structure ?
//	std::vector<size_t> b_tenor(1) ; b_tenor[0] = 1 ;	
//	std::vector<double> strikeATM(1) ; strikeATM[0] =  4.23/100 ; //2.78/100 ;  //1.2341/100 ; 0.36/100 ;
//	std::vector<double> volQuotes(1) ; volQuotes[0] =  16.14/100 ; //29.41/100 ; //70.97/100 ;50.81/100 ;
//	std::vector<double> skew(1) ; skew[0] =  -1.1 ; //-7.4 ; //-3.8 ;
//	double shift = 5.0/10000 ;
//	CourbeInput_PTR courbe_PTR_test(createCourbeInput(7));
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
//	
//	CoTerminalSwaptionVol_CONSTPTR coterminalSwaptionVol_PTR 
//			(new CoTerminalSwaptionVol(	volQuotes, a_expiry,b_tenor, strikeATM)) ;
//	CheyetteDD_CostFunctionLevel_PTR cheyetteBaseCostFunctionLevel_PTR(
//				new CheyetteDD_CostFunctionLevel(coterminalSwaptionVol_PTR, indexSwaption, approx)) ;
//
//	CoTerminalSwaptionSkew_CONSTPTR coTerminalSwaptionSkew_PTR
//			(new CoTerminalSwaptionSkew(skew, a_expiry,b_tenor, strikeATM, shift)) ;
//	CheyetteDD_CostFunctionSkew_PTR cheyetteBaseCostFunctionSkew_PTR(
//				new CheyetteDD_CostFunctionSkew(coTerminalSwaptionSkew_PTR, indexSwaption, approx)) ;
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
//		//CheyetteDD_LocalCalibrator(	const QuantLib::Size & maxIterations,       
//		//						const QuantLib::Real & rootEpsilon,        
//		//						const QuantLib::Real & functionEpsilon,    
//		//						Array sigmaInitiate, Array mInitiate,
//		//						Array calibrated_sigma, Array calibrated_m,
//		//						CheyetteDD_CostFunctionLevel_PTR		costFunctionLevel_PTR,
//		//						CheyetteDD_CostFunctionSkew_PTR		costFunctionSkew_PTR)
//
//	CheyetteDD_LocalCalibrator calibrator(	maxIterations,  rootEpsilon,   functionEpsilon,    
//											sigmaInitiate, mInitiate,
//											calibrated_sigma, calibrated_m,
//											cheyetteBaseCostFunctionLevel_PTR,
//											cheyetteBaseCostFunctionSkew_PTR) ; 
//
//	calibrator.solve() ;
//
//}


//void test_CheyetteCalibrationMarketData()
//{
//	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
//	test_calib_Cheyette_OneFile(mkt_file_list[0]);
//}
//
//void retrieveMarketData(const std::string& mkt_data_file, size_t coterminal)
//{
//	std::string folder_name;
//	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
//	folder_name+=base_name_file;
//	LMMPATH::reset_Output_SubFolder(folder_name );
//
//	size_t model_nbYear = 3 ;
//
//	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketData(model_nbYear, mkt_data_file);
//
//}
//void test_calib_Cheyette_OneFile(const std::string& mkt_data_file)
//{
//
////creation de config Cheyette
//	double k = 0.2 ;
//	double sigmaParam = 0.4 ;
//	double mParam = 0.4 ;
//	Piecewiseconst_RR_Function sigma(model_nbYear, sigmaParam) ; 
//	Piecewiseconst_RR_Function m(model_nbYear, mParam) ;
//
//	CheyetteDD_Model::CheyetteDD_Parameter cheyetteDD_Param(k,sigma, m) ;
//
//	CourbeInput_PTR courbeInput_PTR(new 
//		CourbeInput(upperMatrixZC_To_CourbeInput_yield(pLmmSwaptionMarketData->get_LiborQuotes()) )
//								) ;
//
//	CheyetteDD_CalibrationConfig config(cheyetteDD_Param, courbeInput_PTR) ;
//	config.model_nbYear_ = model_nbYear ;
////calibration
//	marketData_Cheyette_LocalCalibration(config, pLmmSwaptionMarketData);
//
//}
//
//
//void marketData_Cheyette_LocalCalibration(	const CheyetteDD_CalibrationConfig& config,
//											LmmSwaptionMarketData_PTR pLmmSwaptionMarketData)
//{
//	size_t nbYear = pLmmSwaptionMarketData->get_nbYear();
//	std::string base_file_name = pLmmSwaptionMarketData->get_MarketDataBaseFileName();
//
//	Tenor tenorfixedleg = Tenor::_1YR ;
//	Tenor tenorfloatleg = Tenor::_6M  ;
//	size_t fixedfloatRatio = tenorfixedleg.ratioTo(tenorfloatleg);
//
//	std::string base_name;
//	base_name = base_file_name+"_CheyetteDD_LocalCalibration" ;
//
//	//create LMM components
//	LMMTenorStructure_PTR pLMMTenorStructure( new LMMTenorStructure(tenorfloatleg,nbYear) );
//
////creation du Cheyette DD Model
//	size_t shiftChoice = config.shiftChoice_ ;
//	CourbeInput_PTR courbeInput = config.courbeInput_PTR_ ;
//
//	CheyetteDD_Model_PTR cheyette_PTR(new CheyetteDD_Model(courbeInput, config.cheyetteDD_Param_, shiftChoice)) ;
//
////creation du ApproxPricer
//	std::cout << "creation d'une fausse swaption pour approx Pricer Cheyette ; à changer" << std::endl ;
//	double strike          = 0.04;
//	LMM::Index  indexStart = 2; 
//	LMM::Index  indexEnd   = 4; 
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_12M;
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 2) );   
//	
//	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
//
//	VanillaSwaption_PTR swaption_PTR(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR approxPricer(
//									new CheyetteDD_VanillaSwaptionApproxPricer(cheyette_PTR, swaption_PTR));	
//	
////creation CostFunction
//	CoTerminalSwaptionVol_CONSTPTR coTerminalSwaptionVol_PTR(
//			new CoTerminalSwaptionVol(upperMatrixToCoTerminalVol(pLmmSwaptionMarketData->get_SwaptionQuotes_ATM(), 
//																config.model_nbYear_)
//									)) ;
//	CoTerminalSwaptionSkew_CONSTPTR coTerminalSwaptionSkew_PTR(
//					new CoTerminalSwaptionSkew(upperMatrixToCoTerminalSkew(pLmmSwaptionMarketData->get_SwaptionQuotes_skew(),
//													config.model_nbYear_,
//													config.strikeBump_))) ;
//
//	size_t indexSwaption = 0 ;
//	CheyetteDD_CostFunctionLevel_PTR cheyetteCostFunctionLevel_PTR(
//									new CheyetteDD_CostFunctionLevel(coTerminalSwaptionVol_PTR, indexSwaption, approxPricer) ) ;
//	CheyetteDD_CostFunctionSkew_PTR cheyetteCostFunctionSkew_PTR(
//									new CheyetteDD_CostFunctionSkew(coTerminalSwaptionSkew_PTR, indexSwaption, approxPricer) ) ;
//
//// Create Calibrator + solve()
//	size_t arraySize = config.cheyetteDD_Param_.sigma_.gety_().size() ;
//	QuantLib::Array arraySigma(vectorToArray(config.cheyetteDD_Param_.sigma_.getx_())) ;
//	QuantLib::Array arrayM(vectorToArray(config.cheyetteDD_Param_.m_.getx_())) ;
//	
//	CheyetteDD_LocalCalibrator cheyetteCalibrator(500, //maxIter
//												 1e-8,   //x_epsilon (root epsilon)
//												 1e-8,   //f_epsilon (function epsilon)   
//												 arraySigma, arrayM, arraySigma, arrayM,
//												 cheyetteCostFunctionLevel_PTR, cheyetteCostFunctionSkew_PTR) ;
//	cheyetteCalibrator.solve() ;
//
////ecriture output
//	std::ostringstream file_result_stream;file_result_stream<<base_name<<"_result.csv";
//	std::string file_calibration_result(file_result_stream.str());
//	//lmmCalibrator.printPlusPlus(file_calibration_result);
//	
//	//cheyetteCalibrator.printPlusPlus(file_calibration_result) ;
//
//	//std::ostringstream file_gDelegate_stream;file_gDelegate_stream<<base_name<<"_gDelegate.csv";
//	//std::string file_gDelegate_vol(file_gDelegate_stream.str() );
//	//pGMatrixMapping->print(file_gDelegate_vol);
//
//	//std::ostringstream file_vol_stream;file_vol_stream<<base_name<<"_vol.csv";
//	//std::string file_calibrated_vol(file_vol_stream.str() );
//	//pNoShifted_HGVolatilityParam->print( file_calibrated_vol );
//
//
//	//{
//	//	std::string common_result_file_name = "calib_result_CheyetteLocal.csv";
//
//	//	std::string full_common_result_file = LMMPATH::get_Root_OutputPath() + common_result_file_name ;
//
//	//	std::ofstream final_result ;
//	//	final_result.open(full_common_result_file.c_str(), std::ios::app);
//
//	//	final_result<<std::endl<<std::endl<< "============= Test At    "<<LMMPATH::get_TimeDateNow()
//	//		<<",,,,,, Error LInf, "<<lmmCalibrator.get_QuoteError_LInf() <<std::endl ;
//	//	final_result<< lmmCalibrator.get_BaseGeneral_Result_Info();
//
//	//	final_result.close();	
//	//}
//}

QuantLib::Array vectorToArray(std::vector<double> v)
{
	QuantLib::Array a(v.size()) ;
	for (size_t i = 0 ; i < v.size() ; ++i)
	{
		a[i] = v[i] ;
	}
	return a ;
}


//std::vector<double> arrayToVector(QuantLib::Array a)
//{
//	std::vector<double> vect ;
//	vect.reserve(a.size()) ;
//	for (size_t i = 0 ; i < a.size() ; ++i)
//	{
//		vect[i] = a[i] ;
//	}
//	return vect ;
//}
//
////CoTerminalSwaptionVol upperMatrixToCoTerminalVol(UpperTriangleVanillaSwaptionQuotes_PTR upperQuotes,
////												  size_t coterminal)
////{
////	// 1st row and column not used! size = nbYear + 1
////	size_t nbLignes		= upperQuotes->size1() ;	//nb year +1
////	size_t nbColonnes	= upperQuotes->size2() ;	//nb year +1
////	assert(nbLignes == nbColonnes) ;
////
////	for (size_t i = 1 ; i < nbLignes-1 ; ++i)
////	{
////		//boost::numeric::ublas::matrix<SwaptionQuote >
////		UpperTriangularVanillaSwaptionQuotes quotes = upperQuotes->get_UpperTriangularVanillaSwaptionQuotes() ;
////		SwaptionQuote qu = quotes(i, coterminal - i) ;			//typedef std::pair<VanillaSwaption, double> SwaptionQuote; 
////		vectQuotes[i] = qu.second ; 
////		vectExpiry[i] = i ;
////		vectTenor[i] = coterminal - i ;
////		UpperTriangularDoubleMatrix strike = upperQuotes->get_UpperTriangularStrike() ;
////		vectStrike[i] = strike(i, coterminal - i) ;
////	}
////
////	return CoTerminalSwaptionVol(vectQuotes, vectExpiry, vectTenor, vectStrike) ;
////}
//
////CoTerminalSwaptionSkew upperMatrixToCoTerminalSkew(UpperTriangleVanillaSwaptionQuotes_PTR upperQuotes,
////												  size_t coterminal, double shift)
////{
////	// 1st row and column not used! size = nbYear + 1
////	size_t nbLignes		= upperQuotes->size1() ;	//nb year +1
////	size_t nbColonnes	= upperQuotes->size2() ;	//nb year +1
////	assert(nbLignes == nbColonnes) ;
////	
////	std::vector<double> vectQuotes(coterminal), vectStrike(coterminal) ;
////	std::vector<size_t>	vectExpiry(coterminal), vectTenor(coterminal) ;
////	vectQuotes[0] = 0 ;
////	vectExpiry[0] = 0 ;
////	vectTenor[0] = 0 ;
////	vectStrike[0] = 0 ;
////	for (size_t i = 1 ; i < nbLignes-1 ; ++i)
////	{
////		//boost::numeric::ublas::matrix<SwaptionQuote >
////		UpperTriangularVanillaSwaptionQuotes quotes = upperQuotes->get_UpperTriangularVanillaSwaptionQuotes() ;
////		SwaptionQuote qu = quotes(i, coterminal - i) ;			//typedef std::pair<VanillaSwaption, double> SwaptionQuote; 
////		vectQuotes[i] = qu.second ; 
////		vectExpiry[i] = i ;
////		vectTenor[i] = coterminal - i ;
////		UpperTriangularDoubleMatrix strike = upperQuotes->get_UpperTriangularStrike() ;
////		vectStrike[i] = strike(i, coterminal - i) ;
////	}
////	return CoTerminalSwaptionSkew(vectQuotes, vectExpiry, vectTenor, vectStrike, shift) ;
//}