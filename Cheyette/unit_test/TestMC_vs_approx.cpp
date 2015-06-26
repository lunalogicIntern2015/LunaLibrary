#include "TestMC_vs_approx.h"

//pour les tests dans le main()
	//std::vector<double> x, m_y, sigma_y ;
	//x.push_back(0.) ; x.push_back(10.) ; x.push_back(20.) ; 
	//double m_param = 1. ;
	//double k = 0.2 ;
	//double sigm = 0.15 ;
	//m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//sigma_y.push_back(sigm) ; sigma_y.push_back(2. * sigm) ;
	//nbSimus[0] = 5 ; nbSimus[1] = 5 ; nbSimus[2] = 5 ; //nbSimus[3] = 100000 ;// nbSimus[4] = 200000 ;

//courbe 0 : courbe plate à 1%
//courbe 1 : courbe test 
//courbe 2 : courbe marché interbancaire du 22-06-15
CourbeInput_PTR createCourbeInput(int curveChoice)
{
	switch (curveChoice)
	{
	case 1:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0.0 ;
		listeMatu.push_back(0) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(1) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(2) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(3) ;	tauxZC.push_back(1./100 + translation) ;  
		listeMatu.push_back(4) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(5) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(10) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(15) ;	tauxZC.push_back(1./100 + translation) ;  
		listeMatu.push_back(20) ;	tauxZC.push_back(1./100 + translation) ;
		listeMatu.push_back(25) ;	tauxZC.push_back(1./100 + translation) ;
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
		return courbe_PTR_test ;
		break ;
		   }
	case 2:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0.0 ;
		listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100 + translation) ; 
		listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100 + translation) ; 
		listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100 + translation) ; 
		listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100 + translation) ;  
		listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100 + translation) ; 
		listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100 + translation) ; 
		listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100 + translation) ; 
		listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100 + translation) ;  
		listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100 + translation) ;
		listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100 + translation) ;
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
		return courbe_PTR_test ;
		break ;
		   }
	case 3:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0.0 ;
		listeMatu.push_back(0) ;			tauxZC.push_back(-0.12/100  + translation) ; 
		listeMatu.push_back(0.002777778) ;	tauxZC.push_back(-0.12/100  + translation) ; 
		listeMatu.push_back(0.019444444) ;	tauxZC.push_back(-0.131/100  + translation) ; 
		listeMatu.push_back(0.083333333) ;	tauxZC.push_back(-0.064/100  + translation) ; 
		listeMatu.push_back(0.166666667) ;	tauxZC.push_back(-0.039/100  + translation) ;  
		listeMatu.push_back(0.25) ;			tauxZC.push_back(-0.014/100  + translation) ; 
		listeMatu.push_back(0.5) ;			tauxZC.push_back(0.05/100  + translation) ; 
		listeMatu.push_back(0.583333333) ;	tauxZC.push_back(0.045/100  + translation) ; 
		listeMatu.push_back(0.666666667) ;	tauxZC.push_back(0.045/100  + translation) ;  
		listeMatu.push_back(0.75) ;			tauxZC.push_back(0.05/100  + translation) ;
		listeMatu.push_back(0.833333333) ;	tauxZC.push_back(0.056/100  + translation) ;
		listeMatu.push_back(0.916666667) ;	tauxZC.push_back(0.064/100  + translation) ; 
		listeMatu.push_back(1) ;			tauxZC.push_back(0.073/100  + translation) ; 
		listeMatu.push_back(1.5) ;			tauxZC.push_back(0.098/100  + translation) ; 
		listeMatu.push_back(2 ) ;			tauxZC.push_back(0.135/100  + translation) ;  
		listeMatu.push_back(3 ) ;			tauxZC.push_back(0.24/100  + translation) ; 
		listeMatu.push_back(4 ) ;			tauxZC.push_back(0.382/100  + translation) ; 
		listeMatu.push_back(5 ) ;			tauxZC.push_back(0.544/100  + translation) ; 
		listeMatu.push_back(6 ) ;			tauxZC.push_back(0.701/100  + translation) ;  
		listeMatu.push_back(7 ) ;			tauxZC.push_back(0.847/100  + translation) ;
		listeMatu.push_back(8 ) ;			tauxZC.push_back(0.981/100  + translation) ;
		listeMatu.push_back(9 ) ;			tauxZC.push_back(1.1/100  + translation) ;  
		listeMatu.push_back(10) ;			tauxZC.push_back(1.198/100  + translation) ; 
		listeMatu.push_back(11) ;			tauxZC.push_back(1.286/100  + translation) ; 
		listeMatu.push_back(12) ;			tauxZC.push_back(1.359/100  + translation) ; 
		listeMatu.push_back(15) ;			tauxZC.push_back(1.516/100  + translation) ;  
		listeMatu.push_back(20) ;			tauxZC.push_back(1.629/100  + translation) ;
		listeMatu.push_back(25) ;			tauxZC.push_back(1.657/100  + translation) ;
		listeMatu.push_back(30) ;			tauxZC.push_back(1.666/100  + translation) ; 
		listeMatu.push_back(35) ;			tauxZC.push_back(1.677/100  + translation) ; 
		listeMatu.push_back(40) ;			tauxZC.push_back(1.675/100  + translation) ;  
		listeMatu.push_back(45) ;			tauxZC.push_back(1.656/100  + translation) ;
		listeMatu.push_back(50) ;			tauxZC.push_back(1.634/100  + translation) ;
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));		   
		return courbe_PTR_test ;
		break ;
		   }
	default :
		throw "courbe non existante" ;
	}
}
	

//param 
//swaption aY bY
//floating leg tenor vs fixed leg tenor (pas vraiment un tenor mais frequence des flux)
//size(x) = size(m_y) + 1 = size(sigma_y) + 1
void test_approx(double strike, size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
				 int curveChoice, int shiftChoice, 
				 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
				 std::vector<size_t> nbSimus)
{	

	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	size_t indexStart = a / tenor ;
	size_t indexEnd = (a+b) / tenor ;
	double fwdProbaT = (double) b ; 

//courbe spot
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//param MC
	std::vector<size_t>			simulationIndex ;
	std::vector<size_t>			discretizationBetweenDates ;
	simulationIndex.push_back(0) ; 	discretizationBetweenDates.push_back(0) ; 
	simulationIndex.push_back(indexStart) ;	discretizationBetweenDates.push_back(indexStart * 100) ; 
	for (size_t i = 1 ; i <= (indexEnd-indexStart) ; ++i)
	{
		simulationIndex.push_back(indexStart + i) ;
		discretizationBetweenDates.push_back(50) ;
	}

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	//MC_CheyetteDD_GenericSwaptionPricer_PTR mc(new MC_CheyetteDD_GenericSwaptionPricer(modele_test_PTR,
	//														rnGenerator,
	//														floatingLegTenor,
	//														fwdProbaT,
	//														simulationIndex,		
	//														discretizationBetweenDates   ) ) ;

	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR,
																						rnGenerator,
																						floatingLegTenor,
																						fwdProbaT,
																						simulationIndex,		
																						discretizationBetweenDates   ) ) ;

//approx
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

	std::cout << "prixSwaption" << std::endl ;
	double approxPrice = approx.prixSwaptionApproxPiterbarg() ;
	double b_barre = approx.get_buffer_b_barre_() ;
	std::cout << "approxPrice : " << approxPrice << std::endl ;
	std::cout << "b barre :     " << b_barre << std::endl ;
	std::cout << "  " << std::endl ;

//MC
	//double nominal = 1 ;
	//LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, a+b+1)) ;
	//GeneticSwap_CONSTPTR genericSwapTest = InstrumentFactory::createVanillaSwap(strike, indexStart, 
	//												indexEnd, floatingLegTenor, fixedLegTenor,
	//												swapStructure, nominal);

	//size_t maturity = indexStart ; //index
	//GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;

	size_t nbMC = nbSimus.size() ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		//std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;  //Generic
		std::vector<double> resMC = mc->price(swaption_PTR_test, nbSimus[i]) ;  //Vanilla
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
	}
	
	double annuityA0 = approx.swapRateDenominator(0., 0.) ;
	double swapRateS0 = approx.swapRate0() ;

	double blackPrice = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																		swapRateS0, strike, a) ;

	//mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
	//						a, b, genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;

	mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
							a, b, swaption_PTR_test, nbSimus, vectPrixMC, vectICinf, vectICsup) ;

}

//param 
//swaption aY bY
//floating leg tenor vs fixed leg tenor (pas vraiment un tenor mais frequence des flux)
//size(x) = size(m_y) + 1 = size(sigma_y) + 1
void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus)
{	

	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	size_t indexStart = a / tenor ;
	size_t indexEnd = (a+b) / tenor ;
	double fwdProbaT = b ;  //size_t vers double ok

//courbe spot
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//param MC
	std::vector<size_t>			simulationIndex ;
	std::vector<size_t>			discretizationBetweenDates ;
	simulationIndex.push_back(0) ; 	discretizationBetweenDates.push_back(0) ; 
	simulationIndex.push_back(indexStart) ;	discretizationBetweenDates.push_back(indexStart * 100) ; 
	for (size_t i = 1 ; i <= (indexEnd-indexStart) ; ++i)
	{
		simulationIndex.push_back(indexStart + i) ;
		discretizationBetweenDates.push_back(50) ;
	}

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	//MC_CheyetteDD_GenericSwaptionPricer_PTR mc(new MC_CheyetteDD_GenericSwaptionPricer(modele_test_PTR,
	//														rnGenerator,
	//														floatingLegTenor,
	//														fwdProbaT,
	//														simulationIndex,		
	//														discretizationBetweenDates   ) ) ;

	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR,
															rnGenerator,
															floatingLegTenor,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates   ) ) ;

//approx
	double strike = 1000. ; //temporaire avant de mettre strike ATM
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

//calcul du strike ATM pour le swap
	double strikeATM = approx.swapRate0() ;
	swap.set_strike(strikeATM) ;
	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approxATM(modele_test_PTR, swaption_PTR_test_ATM);

	std::cout << "prixSwaption" << std::endl ;
	double approxPrice = approxATM.prixSwaptionApproxPiterbarg() ;
	double b_barre = approx.get_buffer_b_barre_() ;
	std::cout << "approxPrice : " << approxPrice << std::endl ;
	std::cout << "b barre :     " << b_barre << std::endl ;
	std::cout << "  " << std::endl ;

//MC
	//double nominal = 1 ;
	//LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, a+b+1)) ;
	//GeneticSwap_CONSTPTR genericSwapTest = InstrumentFactory::createVanillaSwap(strikeATM, indexStart, 
	//												indexEnd, floatingLegTenor, fixedLegTenor,
	//												swapStructure, nominal);

	//size_t maturity = indexStart ; //index
	//GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;

	size_t nbMC = nbSimus.size() ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		//std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;
		std::vector<double> resMC = mc->price(swaption_PTR_test, nbSimus[i]) ;
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
	}
	
	double annuityA0 = approx.swapRateDenominator(0., 0.) ;
	double swapRateS0 = approx.swapRate0() ;

	double blackPrice = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																	swapRateS0, strike, a) ;

	//mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
	//						a, b, genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;

	mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
							a, b, swaption_PTR_test, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
}