#include "TestCalibrator.h"
#include <cassert>



//swaption aYbY, coterminal = a + b
//calibre sur toutes les swaptions coterminales

void testCalib(size_t fileNumber, size_t coterminal)
{
//lecture du fichier et recuperation des market data
	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();

	std::string mkt_data_file = mkt_file_list[fileNumber] ;

	std::string folder_name;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	size_t model_nbYear = 18 ; //definit jusqu'où est lu VCUB

	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketData(model_nbYear, mkt_data_file);
	
	//LmmSwaptionMarketDataFull_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketDataFull(model_nbYear, mkt_data_file) ;


//conversion des market data de matrices à vecteurs coterminaux
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
	Tenor tenorFixed = Tenor::_1YR ;  
 	std::vector<VanillaSwaption_PTR> vectSwaptions = getVectorSwaptions(tenorFloat, tenorFixed,
																	vectExpiry, vectTenor, vectStrike) ;

	double shift = 5.0/10000 ;
	MarketData_PTR mktData(new MarketData(	vectExpiry, vectTenor, vectStrike, vectVol, vectSkew, shift, 
											vectSwaptions, "Black"));
	

//Cheyette DD model
	double k = 0.02 ;
	double sigma = 0.20 ;
	double m = 0.4 ;	

	Piecewiseconst_RR_Function sigmaFunc(coterminal /*- 1*/, sigma) ; 
	Piecewiseconst_RR_Function mFunc(coterminal /*- 1*/, m) ; 

	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigmaFunc, mFunc) ;
	
	int shiftChoice = 1 ;
	
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
	
	courbeInput_PTR->printHorizontal(o) ;
	o << endl ;

//calibrator
	QuantLib::Size maxIterations	= 500 ;
	QuantLib::Real rootEpsilon		= 1./10000 ;
	QuantLib::Real functionEpsilon	= 1./10000 ;

	//recuperation des quotations vol ATM, ATM+5bp, ATM-5bp pour tracer smile
	UpperTriangleVanillaSwaptionQuotes_PTR volATM_ptr =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATM() ;
	UpperTriangleVanillaSwaptionQuotes_PTR volATMmm_ptr =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMmm() ;
	UpperTriangleVanillaSwaptionQuotes_PTR volATMpp_ptr =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMpp() ;

	//minimisation 1D, on passe en parametre sigma[indexSwaption] ou m[indexSwaption]
	for (size_t indexSwaption = 1 ; indexSwaption < vectSwaptions.size() ; ++indexSwaption)
	{
		o << endl ;
		o << "------ Calibration sur la swaption " 
					<< indexSwaption << "Y" << coterminal - indexSwaption << "Y -------" << std::endl ;
		o << std::endl ;

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

		o << "VOL MARKET ; "
			<< volATMmm_ptr->get_UpperTriangularVanillaSwaptionQuotes()(indexSwaption, coterminal - indexSwaption).second	<< ";"
			<< volATM_ptr->get_UpperTriangularVanillaSwaptionQuotes()(indexSwaption, coterminal - indexSwaption).second		<< ";"
			<< volATMpp_ptr->get_UpperTriangularVanillaSwaptionQuotes()(indexSwaption, coterminal - indexSwaption).second	<<std::endl ;
		o << endl ;
	}
	modele_test_PTR->print(o) ;
	
	//compare market smile and model smile for ITM/ATM/OTM strikes
	//generateSmile(model_nbYear, fileNumber, coterminal, o) ;

	std::vector<size_t> nbSimus(7) ;
	nbSimus[0] = 1000 ; nbSimus[1] = 2000 ; nbSimus[2] = 5000 ; nbSimus[3] = 10000 ; nbSimus[4] = 20000 ;
	nbSimus[5] = 50000 ; nbSimus[6] = 100000 ;

	for (size_t a = 1 ; a < coterminal ; ++a)
	{
		test_approx_ATM(a, coterminal - a, tenorFloat, tenorFixed, modele_test_PTR, nbSimus, o) ;	
	}

	o.close() ;
}


void generateSmile(size_t model_nbYear, size_t fileNumber, size_t coterminal, std::ofstream& o)
{
//	//lecture du fichier et recuperation des market data
//	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList() ;
//	std::string mkt_data_file = mkt_file_list[fileNumber] ;
//
////	LmmSwaptionMarketDataFull_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketDataFull(model_nbYear, mkt_data_file) ;
//
//	//recuperation des quotations vol ATM, ATM+5bp, ATM-5bp pour tracer smile
//	UpperTriangleVanillaSwaptionQuotes_PTR volATM_ptr =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATM() ;
//	UpperTriangleVanillaSwaptionQuotes_PTR volATMmm_ptr =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMmm() ;
//	UpperTriangleVanillaSwaptionQuotes_PTR volATMpp_ptr =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMpp() ;
//
//	UpperTriangleVanillaSwaptionQuotes_PTR volATM_p50 =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMp50() ;
//	UpperTriangleVanillaSwaptionQuotes_PTR volATM_p100 =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMp100() ;
//	UpperTriangleVanillaSwaptionQuotes_PTR volATM_p200 =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMp200() ;
//
//	UpperTriangleVanillaSwaptionQuotes_PTR volATM_m50 =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMm50() ;
//	UpperTriangleVanillaSwaptionQuotes_PTR volATM_m100 =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMm100() ;
//	UpperTriangleVanillaSwaptionQuotes_PTR volATM_m200 =  pLmmSwaptionMarketData->get_SwaptionQuotes_ATMm200() ;
//
//	o	<< "swaption ; volATM - 200bp ; volATM - 100bp ; volATM - 50bp ; volATM - 5bp ; volATM ; "
//		<< "volATM + 5 bp ; volATM + 50 bp ; volATM + 100 bp ; volATM + 200 bp " << std::endl ;
//	for (size_t i = 1 ; i < coterminal ; ++i)
//	{
//		o << "vol ATM swaption " << i << "Y" << coterminal - i << "Y ; "  
//			<< volATM_m200->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second	<< ";"
//			<< volATM_m100->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second		<< ";"
//			<< volATM_m50->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second	<< ";"
//			<< volATMmm_ptr->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second	<< ";"
//			<< volATM_ptr->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second		<< ";"
//			<< volATMpp_ptr->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second	<< ";"
//			<< volATM_p50->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second		<< ";"
//			<< volATM_p100->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second	<< ";"
//			<< volATM_p200->get_UpperTriangularVanillaSwaptionQuotes()(i, coterminal - i).second	<<std::endl ;
//	}
}



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

