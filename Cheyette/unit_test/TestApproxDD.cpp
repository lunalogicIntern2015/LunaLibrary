#include <Cheyette\unit_test\TestApproxDD.h>

#include "TestCalibrator.h"


VanillaSwap_PTR createSwapTest()
{
	double strike          = 0.04;
	LMM::Index  indexStart = 0; 
	LMM::Index  indexEnd   = 2; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_12M;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 2) );   
	
	VanillaSwap_PTR pVanillaSwap (
		new VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure));

	pVanillaSwap->show() ;

	return pVanillaSwap ;
}


VanillaSwaption_PTR createSwaptionTest()
{
	double strike          = 0.01;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_12M;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , indexEnd + 1) );   
	
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	VanillaSwaption_PTR vanillaSwaption_PTR_Test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;

	vanillaSwaption_PTR_Test->show() ;

	return vanillaSwaption_PTR_Test ;
}

void testSwap()
{
	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 8; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 10) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	
	swap.show() ;
}


VanillaSwaption_PTR createSwaption(double strike, LMM::Index  indexStart, LMM::Index  indexEnd, 
								   Tenor floatingLegTenorType, Tenor fixedLegTenorType)
{

	LMMTenorStructure_PTR simulationStructure(
				new LMMTenorStructure(floatingLegTenorType, static_cast<int>(indexEnd / floatingLegTenorType.YearFraction())) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	
	VanillaSwaption_PTR vanillaSwaption_PTR_Test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	vanillaSwaption_PTR_Test->show() ;

	return vanillaSwaption_PTR_Test ;
}

//VanillaSwaption_PTR setSwaptionATM_DD(	const CheyetteDD_Model_PTR modele_test_PTR,
//										const Tenor& floatTenor, const Tenor& fixedTenor, 
//										const size_t a, const size_t b)
//{
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//
//	double strike =  - 1000. ; //temporaire avant de mettre strike ATM
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
//	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
//	
//	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//	Cheyette_SwaptionPricer_LinearApprox approx(modele_test_PTR, swaption_PTR_test);
//
////calcul du strike ATM pour le swap
//	double strikeATM = approx.swapRate0() ;
//	swap.set_strike(strikeATM) ;
//	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//
//	return swaption_PTR_test_ATM ;
//}


Cheyette_SwaptionPricer_LinearApprox_PTR createApproxPricer_PTR(size_t xmax, int numCourbe, 
																double k, double sigmaValue, double mValue)
{
	CheyetteDD_Model_PTR	cheyetteDD_Model_Test_PTR(creeCheyetteDD_Modele_PTR(xmax, numCourbe, 
																				k, sigmaValue, mValue)) ;

	Cheyette_SwaptionPricer_LinearApprox_PTR approxPricerTest_PTR(
							new Cheyette_SwaptionPricer_LinearApprox(cheyetteDD_Model_Test_PTR)  ) ; 
	return approxPricerTest_PTR ;
}

Cheyette_SwaptionPricer_LinearApprox_PTR createApproxPricer_PTR(size_t xmax, CourbeInput_PTR pCourbeInput, 
																double k, double sigmaValue, double mValue)
{
	CheyetteDD_Model_PTR	cheyetteDD_Model_Test_PTR(creeCheyetteDD_Modele_PTR(xmax, pCourbeInput, 
																				k, sigmaValue, mValue)) ;

	Cheyette_SwaptionPricer_LinearApprox_PTR approxPricerTest_PTR(
							new Cheyette_SwaptionPricer_LinearApprox(cheyetteDD_Model_Test_PTR)  ) ; 
	return approxPricerTest_PTR ;
}
//Cheyette_SwaptionPricer_LinearApprox_PTR createLinearApproxPricer_PTR(CheyetteModel_PTR cheyetteModel_PTR)
//{
//	Cheyette_SwaptionPricer_LinearApprox_PTR approxPricerTest_PTR(
//							new Cheyette_SwaptionPricer_LinearApprox(cheyetteModel_PTR)  ) ; 
//	return approxPricerTest_PTR ;
//}

//test qualite approximation prix swaption
void testQualiteApprox(size_t xmax, int numCourbe, double k, double sigmaValue, double mValue)
{
	CheyetteDD_Model_PTR pModel = creeCheyetteDD_Modele_PTR(xmax, numCourbe, k, sigmaValue, mValue) ;
	
	Cheyette_SwaptionPricer_LinearApprox pApprox(pModel); 
	
	//swaption 5Y 5Y de strike ATM 
	double strike = 1000 ;
	LMM::Index  indexStart = 10 ; //5Y
	LMM::Index  indexEnd = 20 ;
	Tenor floatingLegTenorType = Tenor::_6M ;
	Tenor fixedLegTenorType = Tenor::_1YR ;
	VanillaSwaption_PTR pSwaption = createSwaption(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType) ;

	pApprox.price(pSwaption) ;
	pSwaption->set_strike(pApprox.get_buffer_s0_()) ; //pour mettre strike ATM

	size_t nbShifts = 13 ;	
	std::vector<double> shifts_bp(nbShifts) ;			
		
	shifts_bp[0] = -250. ; shifts_bp[1] = -200. ; shifts_bp[2] = -150. ; shifts_bp[3] = -100. ; shifts_bp[4] = -50. ;
	shifts_bp[5] = -25. ; shifts_bp[6] = 0. ; shifts_bp[7] = 25. ;
	shifts_bp[8] = 50. ; shifts_bp[9] = 100. ; shifts_bp[10] = 150. ; shifts_bp[11] = 200. ; shifts_bp[12] = 250. ;

	//APPROXIMATION pour multiple strikes
	std::vector<std::vector<double>> resApprox = pApprox.priceMultipleStrikes(pSwaption, shifts_bp) ;

//ecriture dans fichier
	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_output_path() ;
	fileName_s << directory << "qualite_approximation.csv" ; 
	std::string fileName = fileName_s.str();

	std::ofstream o;

	o.open(fileName,  std::ios::out | std::ios::app );
	o	<<	std::endl;

	pModel->print(o) ;
	pSwaption->show() ;

	helpPrinter("strikes", resApprox[1], o) ;
	helpPrinter("approx", resApprox[0], o) ;
	helpPrinter("vol Black approx", resApprox[2], o) ;

	//MC (proba forward QT)
	std::vector<size_t> nbSimus(1) ;
	nbSimus[0] = 150000 ;
		
	double swapRate0 = pApprox.get_buffer_s0_() ;
	double annuity0 = pApprox.swapRateDenominator(0., 0., 0.) ;
	size_t coterminal = static_cast<size_t>((indexStart + indexEnd)/2.)	;

	printSwaptionMultipleStrikesMC(coterminal, floatingLegTenorType, pModel,
									pSwaption, nbSimus, shifts_bp, annuity0, swapRate0, o) ;
		
	o	<<	std::endl;
	o.close() ;
}

void lancementQualiteApprox()
{
	size_t xmax = 10 ;
	int numCourbe = 1 ;
	double k = 0.02 ;
	std::vector<double> vectSigma(1) ; /* vectSigma[0] = 0.2 ;  vectSigma[0] = 0.8 ; */ vectSigma[0] = 1. ; 
	std::vector<double> vectM(1) ; /* vectM[0] = 0. ; vectM[1] = 0.5 ; */ vectM[0] = 1. ;
	for (size_t i_sigma = 0 ; i_sigma < vectSigma.size() ; ++i_sigma)
	{
		for (size_t i_m = 0 ; i_m < vectM.size() ; ++i_m)
		{
			double sigmaValue = vectSigma[i_sigma] ;
			double mValue = vectM[i_m] ;
			testQualiteApprox(xmax, numCourbe, k, sigmaValue, mValue) ;	
		}		
	}	
}

void test_y_barre()
{
	size_t xmax = 20 ;
	int numCourbe = 1 ;  
	double k(0.02) ;
	double sigmaValue = 0.25 ;
	double mValue = 0. ;

//cas m(t) = 0
	Cheyette_SwaptionPricer_LinearApprox_PTR approx = createApproxPricer_PTR(xmax, numCourbe, 
																			k, sigmaValue, mValue) ;

	double t = 1. ;

	double r0 = approx->get_CheyetteModel()->get_courbeInput_PTR()->get_f_0_t(0) ;
	double res_integrale =  pow(0.25 * r0, 2) * (1. - exp(- 2. * k * t)) / (2. * k) ;

	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;

	VanillaSwaption_PTR vanillaSwaption_PTR_Test = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;

	double prix = approx->price(vanillaSwaption_PTR_Test) ;   //permet d'initialiser y_barre

	std::cout << "t  = " << t << std::endl ;
	std::cout << "integrale_main   : " << res_integrale << std::endl ;
	std::cout << "integrale_classe : " << approx->get_buffer_y_bar_()(t) << std::endl ;
	std::cout << " " << std::endl ;

	t = 10. ;
	indexStart = 20;	//car t = 10 et tenor::_6M
	indexEnd = 40 ;		//quelconque
	VanillaSwaption_PTR vanillaSwaption_PTR_Test2 = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;
	double prix2 = approx->price(vanillaSwaption_PTR_Test2) ;   //permet d'initialiser y_barre


	res_integrale =  pow(0.25 * r0, 2) * (1. - exp(- 2. * k * t)) / (2. * k) ;

	std::cout << "t  = " << t << std::endl ;
	std::cout << "integrale_main   : " << res_integrale << std::endl ;
	std::cout << "integrale_classe : " << approx->get_buffer_y_bar_()(t) << std::endl ;

}

void test_ZC_swapRate_Num_Denom()
{
	size_t xmax = 25 ;
	int numCourbe = 2 ;  
	double k = 0.02 ;
	double sigmaValue = 0.25 ;
	double mValue = 0.25 ;

	CheyetteDD_Model_PTR cheyetteDD_Model_PTR_Test = creeCheyetteDD_Modele_PTR(xmax, numCourbe, 
																				k, sigmaValue, mValue) ;
	cheyetteDD_Model_PTR_Test->show() ;

	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;

	VanillaSwaption_PTR vanillaSwaption_PTR_Test = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;
	vanillaSwaption_PTR_Test->show() ;

	Cheyette_SwaptionPricer_LinearApprox_PTR approxPricerTest_PTR(
			new Cheyette_SwaptionPricer_LinearApprox(cheyetteDD_Model_PTR_Test)) ;
	
//t = 0 avec la courbe spot
	std::cout << "derivee_1_classe T = 1Y   " << approxPricerTest_PTR->ZC_1stDerivative_on_xt(0., 1.0, 0., 0.) << std::endl ;
	std::cout << "derivee_1_main   T = 1Y   " << -exp(- 0.85/100) * (1 - exp(-1)) << std::endl ;
	std::cout << "derivee_1_classe T = 2Y   " << approxPricerTest_PTR->ZC_1stDerivative_on_xt(0., 2.0, 0., 0.) << std::endl ;
	std::cout << "derivee_1_main   T = 2Y   " << -exp(- 2 * 0.9/100) * (1 - exp(-2)) << std::endl ;
	std::cout << "   " << std::endl ;
	std::cout << "derivee_2_classe T = 1Y   " << approxPricerTest_PTR->ZC_2ndDerivative_on_xt(0., 1.0, 0., 0.) << std::endl ;
	std::cout << "derivee_2_main   T = 1Y   "  << exp(- 0.85/100) * pow(1 - exp(-1),2) << std::endl ;
	std::cout << "   " << std::endl ;
//swap rate numerator
	std::cout << "swap rate numerator   " << approxPricerTest_PTR->swapRateNumerator(0., 0., 0.) << std::endl ;
	std::cout << "swap rate numerator   " <<  exp(- 1 * 0.85/100) - exp(- 3 * 0.92/100) << std::endl ;
	std::cout << "   " << std::endl ;
//swap rate denominator
	std::cout << "swap rate denominator " << approxPricerTest_PTR->swapRateDenominator(0., 0., 0.) << std::endl ;
	//somme sur les flux fixes : (2Y et 3Y) (delta_fixed = 1)
	std::cout << "swap rate denominator " << 1 * (exp(- 2 * 0.9/100) + exp(- 3 * 0.92/100)) << std::endl ;  

//t = 0.5
	std::cout << " " << std::endl ;
	std::cout << "------------  ZC  ------------------------------------------" << std::endl ;
	std::cout << " " << std::endl ;

	std::cout << "P(0.5, 1, 0.1, 0.1)  " << cheyetteDD_Model_PTR_Test->P(0.5, 1, 0.1, 0.1)  << std::endl ;
	double g = 1 - exp(-0.5) ;
	std::cout << exp(- 0.85/100) / exp(- 0.825/100 * 1/2.) * exp(- 0.1 * g - 1/2. * 0.1 * g * g) ;
	std::cout << " " << std::endl ;

	std::cout << "P(0.5, 1, 1, 1)      " << cheyetteDD_Model_PTR_Test->P(0.5, 1, 1, 1)  << std::endl ;
	std::cout << exp(- 0.85/100) / exp(- 0.825/100 * 1/2.) * exp(- g - 1/2. * g * g) ;
	std::cout << " " << std::endl ;

	std::cout << "P(0.5, 1, 10, 10)    " << cheyetteDD_Model_PTR_Test->P(0.5, 1, 10, 10)  << std::endl ;
	std::cout << exp(- 0.85/100) / exp(- 0.825/100 * 1/2.) * exp(- 10 * g - 1/2. * 10 * g * g) ;
	std::cout << " " << std::endl ;

//swap rate numerator
	//on fixe arbitrairement x_t = 1 pour le test (t = 0.5)
	std::cout << "------------  t = 0.5  ------------" << std::endl ;
	double y_bar_t = approxPricerTest_PTR->get_buffer_y_bar_()(0.5) ;
	std::cout << "swap rate numerator   " << approxPricerTest_PTR->swapRateNumerator(0.5, 1, y_bar_t) << std::endl ;
	std::cout << "swap rate numerator   " <<  cheyetteDD_Model_PTR_Test->P(0.5, 1, 1, y_bar_t) +
												- cheyetteDD_Model_PTR_Test->P(0.5, 3, 1, y_bar_t) << std::endl ;
	std::cout << "   " << std::endl ;
//swap rate denominator
	std::cout << "swap rate denominator " << approxPricerTest_PTR->swapRateDenominator(0.5, 1, y_bar_t) << std::endl ;
	//somme sur les flux fixes : (2Y et 3Y) (delta_fixed = 1)
	std::cout << "swap rate denominator " << 1 * (cheyetteDD_Model_PTR_Test->P(0.5, 2, 1, y_bar_t) 
												+ cheyetteDD_Model_PTR_Test->P(0.5, 3, 1, y_bar_t)) << std::endl ;  

//t = 2
	std::cout << "P(2, 10, 1, 1)    " << cheyetteDD_Model_PTR_Test->P(2, 10, 1, 1)  << std::endl ;
	g = 1 - exp(- 8) ;
	std::cout << exp(- 1.5/100 * 10) / exp(- 0.9/100 * 2) * exp(- 1 * g - 1/2. * 1 * g * g) ;


	std::cout << " " << std::endl ;
	std::cout << "------------  test swap rate  ------------------------------" << std::endl ;
	std::cout << " " << std::endl ;

	std::cout << "buffer T0 " << approxPricerTest_PTR->get_buffer_T0_() << std::endl ;  //1er flux
	std::cout << "buffer TN " << approxPricerTest_PTR->get_buffer_TN_() << std::endl ;  //dernier flux

	// t = 0.5
	double ge1 = 1 - exp(- 0.5) ;

	double PtT0 = exp( - 1 * 0.85/100) / exp( - 0.5 * 0.825/100) * exp(- 0.2 * ge1 - 1/2. * y_bar_t * ge1 * ge1) ;

	//flux fixe 2Y
	double ge3 = 1 - exp(- (2 - 0.5)) ;
	double PtT2Y = exp( - 2 * 0.9/100) / exp( - 0.5 * 0.825/100) * exp(- 0.2 * ge3 - 1/2. * y_bar_t * ge3 * ge3) ;
	//flux fixe 3Y
	double ge4 = 1 - exp(- (3 - 0.5)) ;
	double PtT3Y = exp( - 3 * 0.92/100) / exp( - 0.5 * 0.825/100) * exp(- 0.2 * ge4 - 1/2. * y_bar_t * ge4 * ge4) ;

	std::cout << "swapRate(0.5, 0.2) : " << approxPricerTest_PTR->swapRate(0.5, 0.2, y_bar_t) << std::endl ;
	std::cout << (PtT0 - PtT3Y) / (PtT2Y + PtT3Y) << std::endl ;

}

void test_swapRate_inverse()
{
	size_t xmax = 4 ;
	int numCourbe = 2 ;  
	double k = 0.02 ;
	double sigmaValue = 0.25 ;
	double mValue = 0.2 ;

	CheyetteDD_Model_PTR cheyetteDD_Model_PTR_Test = creeCheyetteDD_Modele_PTR(xmax, numCourbe, 
																				k, sigmaValue, mValue) ;
	cheyetteDD_Model_PTR_Test->show() ;


	double strike          = 0.04;
	LMM::Index  indexStart = 4; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;

	VanillaSwaption_PTR vanillaSwaption_PTR_Test = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;
	vanillaSwaption_PTR_Test->show() ;

	Cheyette_SwaptionPricer_LinearApprox_PTR approxPricerTest_PTR(
		new Cheyette_SwaptionPricer_LinearApprox(cheyetteDD_Model_PTR_Test)) ;

	std::cout << "approxPricerTest_PTR->swapRate(0, 0) : " << approxPricerTest_PTR->swapRate(0., 0., 0.) <<std::endl ;   //OK
	std::cout << "approxPricerTest_PTR->swapRate(0, 0) : " << approxPricerTest_PTR->swapRate(0.1, 0., 0.) <<std::endl ;  
}

void test_fonction_inverse()
{
	size_t xmax = 2 ;
	int numCourbe = 1 ;  
	double k = 1 ;
	double sigmaValue = 0.25 ;
	double mValue = 0.2 ;

	CheyetteDD_Model_PTR modele_test_PTR = creeCheyetteDD_Modele_PTR(xmax, numCourbe, 
																				k, sigmaValue, mValue) ;

	double strike          = 0.04;
	LMM::Index  indexStart = 4; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;

	VanillaSwaption_PTR swaption_PTR_test = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;	
	Cheyette_SwaptionPricer_LinearApprox_PTR approx 
		(new Cheyette_SwaptionPricer_LinearApprox(modele_test_PTR)) ; 

	//std::cout << "-----------  test fonction inverse  -----------------------" << std::endl ;
	//std::cout << " " << std::endl ;
	//std::cout	<< "swapRate(0, 0) : " << approx.swapRate(0., 0., 0.) 
	//			<< ", approx.get_buffer_s0_ : " << approx.get_buffer_s0_() << std::endl ;
	//std::cout << "inverse : " << approx.inverse(0., approx.swapRate(0, 0)) << " vs 0" << std::endl ;
	//std::cout << " " << std::endl ;
	//std::cout << "swapRate(0.5, 2) : " << approx.swapRate(0.5, 2) << std::endl ;
	//std::cout << "inverse : " << approx.inverse(0.5, approx.swapRate(0.5, 2)) << " vs 2" << std::endl ;
	//std::cout << " " << std::endl ;
	//std::cout << "swapRate(0.5, 0.2) : " << approx.swapRate(0.5, 0.2) << std::endl ;
	//std::cout << "inverse : " << approx.inverse(0.5, approx.swapRate(0.5, 0.2)) << " vs 0.2" << std::endl ;
	//std::cout << " " << std::endl ;	
	//std::cout << "swapRate(1, 0.5) : " << approx.swapRate(1, 0.5) << std::endl ;
	//std::cout << "inverse : " << approx.inverse(1, approx.swapRate(1, 0.5)) << " vs 0.5" << std::endl ;
	//std::cout << " " << std::endl ;
	//std::cout << "swapRate(1, 5) : " << approx.swapRate(1, 5) << std::endl ;
	//std::cout << "inverse : " << approx.inverse(1, approx.swapRate(1, 5)) << " vs 5" << std::endl ;

}


double f_num(double x){return pow(x, 4) ;}
double f_num_x(double x){return 4 * pow(x, 3) ;}
double f_num_x2(double x){return 12 * pow(x, 2) ;}
double f_denom(double x){return pow(x, 2) ;}
double f_denom_x(double x){return 2 * x ;}
double f_denom_x2(double x){return 2 ;}

void test_derivatives()
{
	size_t xmax = 2 ;
	int numCourbe = 1 ;  
	double k = 1 ;
	double sigmaValue = 0.5 ;
	double mValue = 0. ;

	CheyetteDD_Model_PTR modele_test_PTR = creeCheyetteDD_Modele_PTR(xmax, numCourbe, 
																				k, sigmaValue, mValue) ;

	double strike          = 0.04;
	LMM::Index  indexStart = 4; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;

	VanillaSwaption_PTR swaption_PTR_test = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;	
	Cheyette_SwaptionPricer_LinearApprox_PTR approx 
		(new Cheyette_SwaptionPricer_LinearApprox(modele_test_PTR)) ; 


	std::cout << " " << std::endl ;
	std::cout << "------------  d sigma_r / dx  --------------" << std::endl ;
	std::cout << " " << std::endl ;

	std::cout << approx->get_CheyetteModel()->localVol_1stDerivative(0., 0.5, 0.)  << std::endl ;
	std::cout << approx->get_CheyetteModel()->localVol_1stDerivative(0.5, 0.5, 0.)  << std::endl ;
	std::cout << approx->get_CheyetteModel()->localVol_1stDerivative(1., 0.5, 0.)  << std::endl ;

	std::cout << " " << std::endl ;
	std::cout << "------------  test 2nd derivative  --------------" << std::endl ;
	std::cout << " " << std::endl ;
 
	double vx = 10. ;

	double n	= f_num(vx) ;
	double n_1	= f_num_x(vx); 
	double n_2	= f_num_x2(vx); 

	double d	= f_denom(vx);
	double d_1	= f_denom_x(vx);
	double d_2	= f_denom_x2(vx);

	double result = (n_2 * d - n * d_2) * d*d - (n_1 * d - n * d_1) * 2 * d * d_1 ;
	result /= (d*d*d*d);
	std::cout << result << " derivee 2nde de x^2" << std::endl ;

	//std::cout << "swapRate(0.5, 0.2) numerateur 1st deriv :  " << approxPricerTest_PTR->swapRateNumerator_1stDerivative(0.5, 0.2) << std::endl ;
	//std::cout << "swapRate(0.5, 0.2) numerateur deriv mano : " << -(1-exp(- 0.5)) * PtT0 + (1 - exp(-3.5)) * PtTN << std::endl ;

	//std::cout << "swapRate(0.5, 0.2) denom 1st deriv :  " << approxPricerTest_PTR->swapRateDenominator_1stDerivative(0.5, 0.2) << std::endl ;
	//std::cout << "swapRate(0.5, 0.2) denom deriv mano : " <<	- (1-exp(- 1.5)) * PtT2Y 
	//															- (1-exp(- 2.5)) * PtT3Y 
	//															- (1-exp(- 3.5)) * PtTN << std::endl ;

	//double n, np ;
	//n = PtT0 - PtTN;
	//np = -(1-exp(- 0.5)) * PtT0 + (1 - exp(-3.5)) * PtTN;
	//double d, dp ;
	//d = PtT2Y + PtT3Y + PtTN ;
	//dp = - (1-exp(- 1.5)) * PtT2Y - (1-exp(- 2.5)) * PtT3Y - (1-exp(- 3.5)) * PtTN ;

	//std::cout << "swapRate(0.5, 0.2) 1st deriv :  " << approxPricerTest_PTR->swapRate_1stDerivative(0.5, 0.2) << std::endl ;
	//std::cout << "swapRate(0.5, 0.2) 1st deriv :  " << (np * d - n * dp) / (d * d) << std::endl ;

	//std::cout << " " << std::endl ;
}


void test_time_average()
{
	size_t xmax = 4 ;
	int numCourbe = 1 ;  
	double k = 0.2 ;
	double sigmaValue = 0.35 ;
	double mValue = 0. ;

	CheyetteDD_Model_PTR modele_test_PTR = creeCheyetteDD_Modele_PTR(xmax, numCourbe, 
																				k, sigmaValue, mValue) ;
	modele_test_PTR->show() ;

	Cheyette_SwaptionPricer_LinearApprox_PTR approx 
		(new Cheyette_SwaptionPricer_LinearApprox(modele_test_PTR)) ; 

	//valide pour m(t) = 0 uniquement
	//double t = 1 ;
	//std::cout	<< "y barre(t) should be " << sigm * sigm * 0.8/100 * 0.8/100 * (1 - exp(- 2 * k * t))/(2 * k) 
	//			<< " , is " << approx.get_buffer_y_bar_t(t)  
	//			<< std::endl ; 

//swaption tres simplifiée
	double strike          = 0.9/100;  //prendre ATM, proche de S0
	LMM::Index  indexStart = 4; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	VanillaSwaption_PTR swaption_PTR_test = createSwaption(strike, indexStart, indexEnd, 
																floatingLegTenorType, fixedLegTenorType) ;	
	swaption_PTR_test->show() ;

	std::cout << "prixSwaption" << std::endl ;
	std::cout << approx->price(swaption_PTR_test) << std::endl ;
	std::cout << "  " << std::endl ;

}


