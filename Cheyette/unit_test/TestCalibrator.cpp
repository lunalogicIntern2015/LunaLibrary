#include "TestCalibrator.h"
#include <cassert>

void helpPrinter(std::string description, std::vector<double> data, std::ofstream& o)
{
	o << description  << " ; " ;
	for (size_t i = 0 ; i < data.size() ; ++i)
	{
		o << data[i] << " ; " ;
	}
	o << std::endl ;
}

CheyetteMarketData_PTR recoverMarketData(const size_t fileNumber, const size_t coterminal)
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

//conversion des market data de matrices à vecteurs coterminaux à l'aide des fonctions de MarketDataConvertor.h
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
	CheyetteMarketData_PTR mktData(new CheyetteMarketData(	courbeInput_PTR,	vectExpiry, vectTenor, vectStrike, 
															vectVol, vectSkew, shift, vectSwaptions, "Black"));
	return mktData ;
}


CheyetteDD_LocalCalibrator_PTR creeCheyetteDD_Calibrator(size_t fileNumber, size_t coterminal) 
{
//parametres du modele DD : valeurs des parametres avant calibration
	double sigma = 0.2 ;
	double m = 0.2 ;
	double k = 0.02 ;

//recuperation des données de marché à partir d'un fichier
	CheyetteMarketData_PTR pMarketData = recoverMarketData(fileNumber, coterminal) ;
	CourbeInput_PTR pCourbeInput = pMarketData->getCourbeInput() ;

//approx pricer
	size_t xmax = coterminal ;
	Cheyette_SwaptionPricer_Approx_PTR cheyetteApprox_PTR = createApproxPricer_PTR(xmax, pCourbeInput, k, sigma, m) ;
	cheyetteApprox_PTR->get_CheyetteModel()->show() ;

//cost functions
	size_t indexSwaption = 1 ;
	Cheyette_CostFunctionLevel_PTR pCheyette_CostFunctionLevel_PTR (
				new Cheyette_CostFunctionLevel(	cheyetteApprox_PTR, pMarketData, indexSwaption)) ;
	Cheyette_CostFunctionSkew_PTR pCheyette_CostFunctionSkew_PTR (
				new Cheyette_CostFunctionSkew(	cheyetteApprox_PTR, pMarketData, indexSwaption)) ;
//calibrator
	QuantLib::Size maxIterations	= 500 ;
	QuantLib::Real rootEpsilon		= 1./10000 ;
	QuantLib::Real functionEpsilon	= 1./10000 ;

	std::vector<double> v_sigma(coterminal /*- 1*/, sigma) ; 
	QuantLib::Array sigmaInitiate = HelperArray::vectorToArray(v_sigma) ; 

	std::vector<double> v_m(coterminal /*- 1*/, m) ; 
	QuantLib::Array mInitiate = HelperArray::vectorToArray(v_m) ;

	QuantLib::Array calibrated_sigma(sigmaInitiate) ;
	QuantLib::Array calibrated_m(mInitiate) ;

	CheyetteDD_LocalCalibrator_PTR calibrator(	new CheyetteDD_LocalCalibrator(maxIterations,  rootEpsilon,   functionEpsilon,    
											sigmaInitiate, mInitiate,
											calibrated_sigma, calibrated_m,
											pCheyette_CostFunctionLevel_PTR,
											pCheyette_CostFunctionSkew_PTR)) ; 

	return calibrator ;
}

void lancementCalibOneFile()
{
	size_t coterminal = 4 ;
//DD
	//size_t fileNumber = 2 ;	
	//CheyetteDD_LocalCalibrator_PTR calibrator = creeCheyetteDD_Calibrator(fileNumber, coterminal) ;
//Quad
	size_t fileNumberMktData = 14 ; //vcun krw //14 ; 
	size_t fileNumberMktData2 = 0 ; 	
	CheyetteQuad_LocalCalibrator_PTR calibrator = creeCheyetteQuad_Calibrator(fileNumberMktData,
																				fileNumberMktData2,
																				coterminal) ;

	calibrator->calibrate() ; 
	calibrator->getCostFunctionLevel_PTR()->getCheyette_ApproxPricer_PTR()->get_CheyetteModel()->show() ;
}

CheyetteQuad_LocalCalibrator_PTR creeCheyetteQuad_Calibrator(size_t fileNumberMktData, 
															 size_t fileNumberMktData_2, 
															 size_t coterminal) 
{
//parametres du modele DD : valeurs des parametres avant calibration
	double a = 0.02 ;
	double b = 0.02 ;
	double c = 0.01 ;
	double k = 0.02 ;

//recuperation des données de marché à partir d'un fichier
	CheyetteMarketData_PTR pMarketData = recoverMarketData(fileNumberMktData, coterminal) ;
	CourbeInput_PTR pCourbeInput = pMarketData->getCourbeInput() ;

//approx pricer
	size_t	xmax = coterminal ;
	//Cheyette_SwaptionPricer_Approx_PTR cheyetteApprox_PTR = createApproxPricer_PTR(xmax, pCourbeInput, k, sigma, m) ;
	Cheyette_SwaptionPricer_QuadApprox_PTR cheyetteApprox_PTR = createQuadApproxPricer_PTR(xmax, pCourbeInput, 
																					k, a, b, c) ;
	cheyetteApprox_PTR->get_CheyetteModel()->show() ;

//cost functions
	size_t indexSwaption = 1 ;
	Cheyette_CostFunctionLevel_PTR pCheyette_CostFunctionLevel_PTR (
				new Cheyette_CostFunctionLevel(	cheyetteApprox_PTR, pMarketData, indexSwaption)) ;
	Cheyette_CostFunctionSkew_PTR pCheyette_CostFunctionSkew_PTR (
				new Cheyette_CostFunctionSkew(	cheyetteApprox_PTR, pMarketData, indexSwaption)) ;

	//cheyette market data 2
	CheyetteMarketData_2_PTR cheyetteMarketData_2_PTR(new CheyetteMarketData_2()) ;
	cheyetteMarketData_2_PTR->init(fileNumberMktData_2) ;
	Cheyette_CostFunctionConvexity_PTR pCheyette_CostFunctionConvexity_PTR (
				new Cheyette_CostFunctionConvexity(	cheyetteApprox_PTR, pMarketData, indexSwaption, 
													cheyetteMarketData_2_PTR)) ;

//calibrator
	QuantLib::Size maxIterations	= 500 ;
	QuantLib::Real rootEpsilon		= 1./10000 ;
	QuantLib::Real functionEpsilon	= 1./10000 ;

	std::vector<double> v_a(coterminal /*- 1*/, a) ; 
	QuantLib::Array aInitiate = HelperArray::vectorToArray(v_a) ; 

	std::vector<double> v_b(coterminal /*- 1*/, b) ; 
	QuantLib::Array bInitiate = HelperArray::vectorToArray(v_b) ;

	std::vector<double> v_c(coterminal /*- 1*/, c) ; 
	QuantLib::Array cInitiate = HelperArray::vectorToArray(v_c) ;

	QuantLib::Array calibrated_a(aInitiate) ;
	QuantLib::Array calibrated_b(bInitiate) ;
	QuantLib::Array calibrated_c(cInitiate) ;

	CheyetteQuad_LocalCalibrator_PTR calibrator(	new CheyetteQuad_LocalCalibrator(maxIterations,  rootEpsilon,   functionEpsilon,    
											aInitiate, bInitiate, cInitiate,
											calibrated_a, calibrated_b, calibrated_c,
											pCheyette_CostFunctionLevel_PTR,
											pCheyette_CostFunctionSkew_PTR,
											pCheyette_CostFunctionConvexity_PTR)) ; 

	return calibrator ;
}


//calibre sur toutes les swaptions coterminales
CheyetteModel_PTR getCalibratedModel(size_t fileNumber, size_t coterminal)
{
//Displaced Diffusion
	CheyetteDD_LocalCalibrator_PTR calibrator = creeCheyetteDD_Calibrator(fileNumber, coterminal) ;

//Quadratic volatility
	//size_t fileNumberMktData_2 = 0 ;
	//CheyetteQuad_LocalCalibrator_PTR calibrator = creeCheyetteQuad_Calibrator(fileNumber, fileNumberMktData_2, coterminal) ;
	
	calibrator->calibrate() ;
	CheyetteModel_PTR calibratedModel = calibrator->getCostFunctionLevel_PTR()->getCheyette_ApproxPricer_PTR()->get_CheyetteModel() ;
	
	calibratedModel->show() ;
	return calibratedModel ;	
}

//swaption aYbY, coterminal = a + b
//calibre sur toutes les swaptions coterminales

void printAllResults_calibratedData(size_t fileNumber, size_t coterminal)
{
	CheyetteMarketData_PTR marketData = recoverMarketData(fileNumber, coterminal) ;
	Tenor floatTenor = Tenor::_6M ;

//ecriture dans fichier
	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_output_path() ;
	fileName_s << directory << "10Y_test_RESULTATS_POUR_RAPPORT.csv" ; 
	std::string fileName = fileName_s.str();

	ofstream o;

	o.open(fileName,  ios::out | ios::app );
	o	<<	endl;
	
//calibrator
	CheyetteModel_PTR calibratedModel = getCalibratedModel(fileNumber, coterminal) ;
	
	calibratedModel->print(o) ;

	std::vector<VanillaSwaption_PTR> vectSwaptions = marketData->getVect_swaptions() ;

	//approx
	Cheyette_SwaptionPricer_Approx_PTR pApproxPricer(new Cheyette_SwaptionPricer_LinearApprox(calibratedModel) ) ;
	//Cheyette_SwaptionPricer_Approx_PTR pApproxPricer(new Cheyette_SwaptionPricer_QuadApprox(calibratedModel) ) ;

	//swaption 0 non definie
	for (size_t a = 10 ; a < coterminal ; ++a)
	{
		////pour les swaptions 1Y, 3Y, 5Y, 10Y, 15Y seulement, commenter sinon

		//bool boolean_a = a == 1 || a == 3 || a == 5 || a == 10 || a == 15 ;
		//size_t b = coterminal - a ;
		//bool boolean_cot =	b == 1 || b == 3 || b == 5 || b == 10 || b == 15 ;
		//bool boolean = boolean_a && boolean_cot ;
		//if (!boolean) continue ;

		o << endl ;
		o   << "SWAPTION " << a << "Y ," << coterminal - a << "Y" << std::endl;
		o   << "startdate ;" << vectSwaptions[a]->getUnderlyingSwap().get_StartDate() << std::endl;
		o   << "enddate ;" << vectSwaptions[a]->getUnderlyingSwap().get_EndDate() << std::endl;
		double strike_0 = vectSwaptions[a]->get_strike() ;		
		o   << "strike ATM ;" << strike_0 << std::endl;
		double approxPrice = pApproxPricer->price(vectSwaptions[a]) ;
		o   << "prix par approximation ;" << approxPrice << std::endl;


		//liste de shifts EN BP à appliquer sur le strike
		size_t nbShifts = 13 ;	
		std::vector<double> shifts_bp(nbShifts) ;			
		
		//shifts_bp[0] = -250. ; shifts_bp[1] = -200. ; 
		//shifts_bp[0] = -150. ; shifts_bp[1] = -100. ; shifts_bp[2] = -50. ;
		//shifts_bp[3] = -25. ; shifts_bp[4] = 0. ; shifts_bp[5] = 25. ;
		//shifts_bp[6] = 50. ; shifts_bp[7] = 100. ; shifts_bp[8] = 150. ; shifts_bp[9] = 200. ; shifts_bp[10] = 250. ;

		shifts_bp[0] = -250. ; shifts_bp[1] = -200. ; shifts_bp[2] = -150. ; shifts_bp[3] = -100. ; shifts_bp[4] = -50. ;
		shifts_bp[5] = -25. ; shifts_bp[6] = 0. ; shifts_bp[7] = 25. ;
		shifts_bp[8] = 50. ; shifts_bp[9] = 100. ; shifts_bp[10] = 150. ; shifts_bp[11] = 200. ; shifts_bp[12] = 250. ;

		//APPROXIMATION pour multiple strikes
		std::vector<std::vector<double>> resApprox = pApproxPricer->priceMultipleStrikes(vectSwaptions[a], shifts_bp) ;
	
		helpPrinter("strikes", resApprox[1], o) ;
		helpPrinter("approx", resApprox[0], o) ;
		helpPrinter("vol Black approx", resApprox[2], o) ;

		//MC (proba forward QT)
		std::vector<size_t> nbSimus(1) ;
		nbSimus[0] = 60000 ;
		
		double swapRate0 = pApproxPricer->get_buffer_s0_() ;
		double annuity0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
		
		printSwaptionMultipleStrikesMC(coterminal, floatTenor, calibratedModel,
										vectSwaptions[a], nbSimus, shifts_bp, annuity0, swapRate0, o) ;

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

		//MC annuity2 (linearise)

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


void printSwaptionMultipleStrikesMC(size_t coterminal, Tenor floatTenor, 
									CheyetteModel_PTR calibratedModel,
									VanillaSwaption_PTR pSwaption, 
									std::vector<size_t> nbSimus, 
									std::vector<double> shifts_bp, double annuity0, double swapRate0,
									std::ofstream& o)
{
		o	<<	endl;
		o   << "Monte Carlo sous proba forward Q T ; T = ;" << coterminal << std::endl;
		o << "nb simus MC ; " << nbSimus[0] << std::endl ;

		LMMTenorStructure_PTR	pTenorStructure(new LMMTenorStructure(floatTenor, coterminal + 1)) ;
		size_t					fwdProbaT = coterminal ;
		size_t					discretizationBetweenDates = 200 ;
		MC_Cheyette_VanillaSwaptionPricer_PTR mcPricer = 
				creeMC_SwaptionPricer_PTR(calibratedModel, pTenorStructure, fwdProbaT, discretizationBetweenDates) ;

		double T0 = pSwaption->getUnderlyingSwap().get_StartDate() ;
		std::cout << T0 << std::endl ;

		for (size_t indexSimu = 0 ; indexSimu < nbSimus.size() ; ++indexSimu)
		{
			std::vector<std::vector<double>> resMC = mcPricer->priceMultipleStrikes(pSwaption, nbSimus[indexSimu], shifts_bp) ;
			helpPrinter("prixSwaptionsPlusieursStrikes", resMC[0], o) ;
			helpPrinter("IC_left", resMC[1], o) ;
			helpPrinter("IC_right", resMC[2], o) ;
			helpPrinter("strikes", resMC[3], o) ;

			std::vector<double> volImp(resMC[0].size()) ;
			for (size_t i = 0 ; i < resMC[0].size() ; i++)
			{
				double vol = NumericalMethods::Black_impliedVolatility(resMC[0][i] / annuity0, swapRate0, resMC[3][i], T0) ;
				volImp[i] = vol ;
			}
			helpPrinter("vol implicite", volImp, o) ;
		}
}

//void lancementAuto()
//{
//
//	
//	std::vector<size_t> vectCoterminal(13) ;
//	vectCoterminal[0] = 2 ; vectCoterminal[1] = 4 ; vectCoterminal[2] = 6 ; vectCoterminal[3] = 8 ;
//	vectCoterminal[4] = 10 ;	vectCoterminal[5] = 11 ; vectCoterminal[6] = 13 ; vectCoterminal[7] = 15 ;
//	vectCoterminal[8] = 16 ; 	vectCoterminal[9] = 18 ; vectCoterminal[10] = 20 ; vectCoterminal[11] = 25 ;
//	vectCoterminal[12] = 30 ;
//
//	//size_t fileNumber = 1 ;
//	for (size_t fileNumber = 12 ; fileNumber < 13 ; ++fileNumber)
//	{
//		for (size_t i = 0 ; i < vectCoterminal.size() ; ++i)
//		{
//			printAllResults_calibratedData(fileNumber, vectCoterminal[i]) ;
//		}
//	}
//}

	//file_list.push_back("VCUB_2008_01_02.csv");// 0
	//	file_list.push_back("VCUB_2008_04_02.csv");// 1  
	//	file_list.push_back("VCUB_2008_07_02.csv");// 2
	//	file_list.push_back("VCUB_2008_10_01.csv");// 3
	//	file_list.push_back("VCUB_2009_01_07.csv");// 4
	//	file_list.push_back("VCUB_2009_04_01.csv");// 5
	//	file_list.push_back("VCUB_2009_07_01.csv");// 6
	//	file_list.push_back("VCUB_2009_10_07.csv");// 7
	//	file_list.push_back("VCUB_2010_01_06.csv");// 8
	//	file_list.push_back("VCUB_2010_04_07.csv");// 9
	//	file_list.push_back("VCUB_2010_07_07.csv");// 10
	//	file_list.push_back("VCUB_2010_10_06.csv");// 11
	//	file_list.push_back("VCUB_2011_01_05.csv");// 12
	//	file_list.push_back("VCUB_2011_04_06.csv");// 13
	//	file_list.push_back("VCUB_2011_07_06.csv");// 14
	//	file_list.push_back("VCUB_2011_10_05.csv");// 15
	//	file_list.push_back("VCUB_2012_01_04.csv");// 16
	//	file_list.push_back("VCUB_2012_04_04.csv");// 17
	//	file_list.push_back("VCUB_2012_07_04.csv");// 18
	//	file_list.push_back("VCUB_2012_10_03.csv");// 19
	//	file_list.push_back("VCUB_2013_01_02.csv");// 20
	//	file_list.push_back("VCUB_2013_04_04.csv");// 21
	//	file_list.push_back("VCUB_2013_07_03.csv");// 22
	//	file_list.push_back("VCUB_2013_10_03.csv");// 23
	//	file_list.push_back("VCUB_2014_01_01.csv");// 24
	//	file_list.push_back("VCUB_2014_02_05.csv");// 25
	//	file_list.push_back("VCUB_2014_03_05.csv");// 26
	//	file_list.push_back("VCUB_2014_04_02.csv");// 27
	//	file_list.push_back("VCUB_2014_05_01.csv");// 28
	//	file_list.push_back("VCUB_2014_06_04.csv");// 29
	//	file_list.push_back("VCUB_2014_07_02.csv");// 30
	//	file_list.push_back("VCUB_2014_08_01.csv");// 31
	//	file_list.push_back("VCUB_2014_08_04.csv");// 31
	//	file_list.push_back("VCUB_2014_08_05.csv");// 33
	//	file_list.push_back("VCUB_2014_08_06.csv");// 34
	//	file_list.push_back("VCUB_2014_08_07.csv");// 35
	//	file_list.push_back("VCUB_2014_08_08.csv");// 36
	//	file_list.push_back("VCUB_2014_08_14.csv");// 37





