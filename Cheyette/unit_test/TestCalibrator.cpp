#include "TestCalibrator.h"
#include <cassert>



//swaption aYbY, coterminal = a + b
//calibre sur toutes les swaptions coterminales

void printAllResults_calibratedData(size_t fileNumber, size_t coterminal)
{
//lecture du fichier et recuperation des market data
	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();

	std::string mkt_data_file = mkt_file_list[fileNumber] ;

	std::string folder_name;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	size_t model_nbYear = coterminal + 2 ; //definit jusqu'où est lu VCUB

	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketData(model_nbYear, mkt_data_file);

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
	//parametres initiaux du modele à calibrer
	double k = 0.02 ;
	double sigma = 1.0 ;
	double m = 0.10 ;	

	Piecewiseconst_RR_Function sigmaFunc(coterminal /*- 1*/, sigma) ;	//modifier pour mettre 15Y 15Y
	Piecewiseconst_RR_Function mFunc(coterminal /*- 1*/, m) ; 

	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigmaFunc, mFunc) ;
	
	int shiftChoice = 1 ;
	
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbeInput_PTR, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//ecriture dans fichier
	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_output_path() ;
	fileName_s << directory << "results_ALL.csv" ; 
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
		QuantLib::Array sigmaInitiate = vectorToArray(v_sigma) ; 

		std::vector<double> v_m(1) ; v_m[0] = m ;
		QuantLib::Array mInitiate = vectorToArray(v_m) ;

		QuantLib::Array calibrated_sigma(sigmaInitiate) ;
		QuantLib::Array calibrated_m(mInitiate) ;

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

	//for (size_t a = 1 ; a < coterminal ; ++a)
	//{
	//	test_approx_ATM(a, coterminal - a, tenorFloat, tenorFixed, modele_test_PTR, nbSimus, o) ;	
	//}

	for (size_t a = 1 ; a < coterminal ; ++a)
	{
		//pour les swaptions 1Y, 3Y, 5Y, 10Y, 15Y seulement, commenter sinon
		bool boolean_a = a == 1 || a == 3 || a == 5 || a == 10 || a == 15 ;
		size_t b = coterminal - a ;
		bool boolean_cot =	b == 1 || b == 3 || b == 5 || b == 10 || b == 15 ;
		bool boolean = boolean_a && boolean_cot ;
		if (!boolean) continue ;

		//swaption 0 non definie
		CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer(
			new CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, vectSwaptions[a])) ;
		o << endl ;
		o   << "SWAPTION " << a << "Y ," << coterminal - a << "Y" << std::endl;
		o   << "startdate ;" << pApproxPricer->get_buffer_UnderlyingSwap_().get_StartDate() << std::endl;
		o   << "enddate ;" << pApproxPricer->get_buffer_UnderlyingSwap_().get_EndDate() << std::endl;
		
		//approx
		double approxPrice = pApproxPricer->prixSwaptionApproxPiterbarg() ;
		o   << "prix par approximation ;" << approxPrice << std::endl;

		//MC forward QT
		double strikeATM_Bloomberg = vectSwaptions[a]->get_strike() ;
		double annuity0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
		double swapRate0 = pApproxPricer->get_buffer_s0_() ;
 
		double sigma_ATM = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuity0, 															
							swapRate0, strikeATM_Bloomberg, a) ;
		//std::vector<size_t> nbSimus(16) ;
		//nbSimus[0] = 1000 ; nbSimus[1] = 3000 ; nbSimus[2] = 5000 ; nbSimus[3] = 10000 ; nbSimus[4] = 30000 ;
		//nbSimus[5] = 50000 ; nbSimus[6] = 100000 ; nbSimus[7] = 300000 ; nbSimus[8] = 500000 ; nbSimus[9] = 800000 ;
		//nbSimus[10] = 1000000 ; nbSimus[11] = 1200000 ; nbSimus[12] = 1500000 ; nbSimus[13] = 1700000 ; nbSimus[14] = 2000000 ;
		//nbSimus[15] = 2500000 ;

		std::vector<size_t> nbSimus(1) ;
		nbSimus[0] = 60000 ;
		o	<<	endl;
		o   << "Monte Carlo sous proba forward Q T ; T = ;" << coterminal << std::endl;
		//o   << "nbSimus ; MC swaption ; IC_left ; IC_right" << std::endl;
		o << "nb simus MC ; " << nbSimus[0] << std::endl ;

		for (size_t indexSimu = 0 ; indexSimu < nbSimus.size() ; ++indexSimu)
		{
			//print_MCforward(nbSimus[indexSimu], pApproxPricer, o) ;
			//std::vector<double> resMC = MCforward(nbSimus[indexSimu], pApproxPricer) ;
			//o << nbSimus[indexSimu] << " ; " << resMC[0] << " ; " << resMC[1] << " ; " << resMC[2] << std::endl ;

			//res[0] = prixSwaptionsPlusieursStrikes ;
			//res[1] = IC_left ; //res[2] = IC_right ; 
			//res[3] = strikes ; //res[4] = moneyness ;

			std::vector<std::vector<double>> resMC = MCforward_MultipleStrikes(nbSimus[indexSimu], pApproxPricer, sigma_ATM) ;
			helpPrinter("prixSwaptionsPlusieursStrikes", resMC[0], o) ;
			helpPrinter("IC_left", resMC[1], o) ;
			helpPrinter("IC_right", resMC[2], o) ;
			helpPrinter("strikes", resMC[3], o) ;
			helpPrinter("moneyness", resMC[4], o) ;

			std::vector<double> volImp(resMC[0].size()) ;
			for (size_t i = 0 ; i < resMC[0].size() ; i++)
			{
					//double vol = NumericalMethods::Black_SwaptionImpliedVolatility(resMC[0][i], annuity0, 															
					//		swapRate0, resMC[3][i], a) ;
					double vol = NumericalMethods::Black_impliedVolatility(resMC[0][i] / annuity0, 															
							swapRate0, resMC[3][i], a) ;
					volImp[i] = vol ;
			}
			helpPrinter("vol implicite", volImp, o) ;

			//strike varie
			//print_smile(a, coterminal - a, strikeATM_Bloomberg, resMC[0], annuity0, swapRate0, o) ;
		}

		//MC annuity
		//o	<<	endl;
		//o   << "Monte Carlo sous proba annuity" << std::endl;
		//o   << "nbSimus ; MC swaption ; IC_left ; IC_right" << std::endl;
		//for (size_t indexSimu = 0 ; indexSimu < nbSimus.size() ; ++indexSimu)
		//{
		//	//print_MCannuity(nbSimus[indexSimu], pApproxPricer, o) ;
		//	std::vector<double> resMC = MCannuity(nbSimus[indexSimu], pApproxPricer) ;
		//	o << nbSimus[indexSimu] << " ; " << resMC[0] << " ; " << resMC[1] << " ; " << resMC[2] << std::endl ;

		//	print_smile(a, coterminal - a, strikeATM_Bloomberg, resMC[0], annuity0, swapRate0, o) ;
		//}

		////MC annuity2 (libnearise)
		//o	<<	endl;
		//o   << "Monte Carlo sous proba annuity (eq 2 lineaire)" << std::endl;
		//o   << "nbSimus ; MC swaption ; IC_left ; IC_right" << std::endl;
		//for (size_t indexSimu = 0 ; indexSimu < nbSimus.size() ; ++indexSimu)
		//{
		//	//print_MCannuity2(nbSimus[indexSimu], pApproxPricer, o) ;
		//	std::vector<double> resMC = MCannuity2(nbSimus[indexSimu], pApproxPricer) ;
		//	o << nbSimus[indexSimu] << " ; " << resMC[0] << " ; " << resMC[1] << " ; " << resMC[2] << std::endl ;

		//	print_smile(a, coterminal - a, strikeATM_Bloomberg, resMC[0], annuity0, swapRate0, o) ;
		//}
		//o	<<	endl;
	}

	o.close() ;

}

void lancementAuto()
{
	std::vector<size_t> vectCoterminal(13) ;
	vectCoterminal[0] = 2 ; vectCoterminal[1] = 4 ; vectCoterminal[2] = 6 ; vectCoterminal[3] = 8 ;
	vectCoterminal[4] = 10 ;	vectCoterminal[5] = 11 ; vectCoterminal[6] = 13 ; vectCoterminal[7] = 15 ;
	vectCoterminal[8] = 16 ; 	vectCoterminal[9] = 18 ; vectCoterminal[10] = 20 ; vectCoterminal[11] = 25 ;
	vectCoterminal[12] = 30 ;

	//size_t fileNumber = 1 ;
	for (size_t fileNumber = 5 ; fileNumber < 30 ; ++fileNumber)
	{
		for (size_t i = 0 ; i < vectCoterminal.size() ; ++i)
		{
			printAllResults_calibratedData(fileNumber, vectCoterminal[i]) ;
		}
	}
}


//calibre sur toutes les swaptions coterminales
CheyetteDD_Model_PTR getCalibratedModel(size_t fileNumber, size_t coterminal)
{
//initial parameters of Cheyette DD model
	double k = 0.02 ;
	double sigma = 0.2 ; //1.0 ;
	double m = 0.4 ; // 0.10 ;	

	Piecewiseconst_RR_Function sigmaFunc(coterminal /*- 1*/, sigma) ; 
	Piecewiseconst_RR_Function mFunc(coterminal /*- 1*/, m) ; 

	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigmaFunc, mFunc) ;
	
	int shiftChoice = 1 ;

	
//lecture du fichier et recuperation des market data
	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();

	std::string mkt_data_file = mkt_file_list[fileNumber] ;

	std::string folder_name;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	size_t model_nbYear = 18 ; //definit jusqu'où est lu VCUB

	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData = get_LmmSwaptionMarketData(model_nbYear, mkt_data_file);


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
		QuantLib::Array sigmaInitiate = vectorToArray(v_sigma) ; 

		std::vector<double> v_m(1) ; v_m[0] = m ;
		QuantLib::Array mInitiate = vectorToArray(v_m) ;

		QuantLib::Array calibrated_sigma(sigmaInitiate) ;
		QuantLib::Array calibrated_m(mInitiate) ;

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
	
	//print le modele calibre
	return modele_test_PTR ;
	
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

