#include <Cheyette\unit_test\TestApproxDD.h>



//CourbeInput_PTR createCourbeInput()
//{
//	std::vector<double> listeMatu, tauxZC ;
//	double translation = 0.0 ;
////courbe des taux plate
//	std::cout << "courbe des taux plate à 1 %" << std::endl ;
//	listeMatu.push_back(0) ;	tauxZC.push_back(1./100 + translation) ; 
//	listeMatu.push_back(1) ;	tauxZC.push_back(1./100 + translation) ; 
//	listeMatu.push_back(2) ;	tauxZC.push_back(1./100 + translation) ; 
//	listeMatu.push_back(3) ;	tauxZC.push_back(1./100 + translation) ;  
//	listeMatu.push_back(4) ;	tauxZC.push_back(1./100 + translation) ; 
//	listeMatu.push_back(5) ;	tauxZC.push_back(1./100 + translation) ; 
//	listeMatu.push_back(10) ;	tauxZC.push_back(1./100 + translation) ; 
//	listeMatu.push_back(15) ;	tauxZC.push_back(1./100 + translation) ;  
//	listeMatu.push_back(20) ;	tauxZC.push_back(1./100 + translation) ;
//	listeMatu.push_back(25) ;	tauxZC.push_back(1./100 + translation) ;
//	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
//
//	//listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100 + translation) ; 
//	//listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100 + translation) ; 
//	//listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100 + translation) ; 
//	//listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100 + translation) ;  
//	//listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100 + translation) ; 
//	//listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100 + translation) ; 
//	//listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100 + translation) ; 
//	//listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100 + translation) ;  
//	//listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100 + translation) ;
//	//listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100 + translation) ;
//	//CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
//
//	//std::cout << "interpolation des taux ZC : " << std::endl ;							//OK
//	//std::cout << courbe_PTR_test->get_tauxZC0(0) << "vs " <<  0.8/100 << std::endl ;
//	//std::cout << courbe_PTR_test->get_tauxZC0(2.5) << "vs " <<  0.91/100 << std::endl ;
//	//std::cout << courbe_PTR_test->get_tauxZC0(20) << "vs " <<  2.5/100 << std::endl ;
//
//	//std::cout << courbe_PTR_test->get_tauxZC0(-1) << "erreur" << std::endl ;
//	//std::cout << courbe_PTR_test->get_tauxZC0(30) << "erreur" << std::endl ;
//
//	//std::vector<double> listeMatu, tauxZC ;
//	//double translation = 0.0 ;
//	//listeMatu.push_back(0) ;	tauxZC.push_back(4.2/100 + translation) ; 
//	//listeMatu.push_back(1) ;	tauxZC.push_back(4.2/100 + translation) ; 
//	//listeMatu.push_back(5) ;	tauxZC.push_back(5.3/100 + translation) ; 
//	//listeMatu.push_back(10) ;	tauxZC.push_back(6.0/100 + translation) ;  
//	//CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
//
//	//std::cout << "taux forward instantane f(0,t) : " << std::endl ;						//OK
//	//std::cout << courbe_PTR_test->get_f_0_t(0.5)  << std::endl ;
//	//std::cout << courbe_PTR_test->get_f_0_t(3.5)  << std::endl ;
//	//std::cout << courbe_PTR_test->get_f_0_t(7)    << std::endl ;
//
//	return courbe_PTR_test ;
//}

CheyetteDD_Model_PTR createCheyetteDD_Model()
{
	int curveChoice = 1 ;
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));

	std::vector<double> x, sigma_y, m_y ;
	x.push_back(0) ; x.push_back(1) ; x.push_back(2) ; 
	m_y.push_back(0.25) ; m_y.push_back(0.5) ;
	sigma_y.push_back(0.25) ; sigma_y.push_back(0.5) ;

	Piecewiseconst_RR_Function sigma = Piecewiseconst_RR_Function(x, sigma_y) ; 
	Piecewiseconst_RR_Function m = Piecewiseconst_RR_Function(x, m_y) ; 
	double k = 0.5 ;

	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR cheyetteDD_Model_PTR_Test(new CheyetteDD_Model( courbe_PTR_test, 
														CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m), shiftChoice )) ;
	cheyetteDD_Model_PTR_Test->show() ;
	return cheyetteDD_Model_PTR_Test ;
}

VanillaSwaption_PTR createSwap()
{
	double strike          = 0.04;
	LMM::Index  indexStart = 0; 
	LMM::Index  indexEnd   = 2; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_12M;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 2) );   
	
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);

	std::cout << "get_EndDate " << simulationStructure->get_EndDate() << std::endl ;
	std::cout << "get_horizon " << simulationStructure->get_horizon() << std::endl ;
	std::cout << "get_nbLIBOR " << simulationStructure->get_nbLIBOR() << std::endl ;
	std::cout << "get_tenorType " << simulationStructure->get_tenorType() << std::endl ;

	std::vector<LMM::Index> vect_fixedLegPaymentIndexSchedule, vect_floatLegPaymentIndexSchedule ;
	vect_fixedLegPaymentIndexSchedule = swap.get_fixedLegPaymentIndexSchedule();
	vect_floatLegPaymentIndexSchedule = swap.get_floatingLegPaymentIndexSchedule() ;
	std::cout << "indices des flux fixes" << std::endl ;
	for (size_t i= 0 ; i < vect_fixedLegPaymentIndexSchedule.size() ; ++i)
	{
		std::cout << vect_fixedLegPaymentIndexSchedule[i] << std::endl ;
	}		//retourne 2 4 6 8
	std::cout << "indices des flux flottants" << std::endl ;
	for (size_t i= 0 ; i < vect_floatLegPaymentIndexSchedule.size() ; ++i)
	{
		std::cout << vect_floatLegPaymentIndexSchedule[i] << std::endl ;
	}		//retourne 1 2 3 4 5 6 7 8

	std::vector<double> vect_tenor_Dates = simulationStructure->get_tenorDate() ;
	for (size_t i = 0 ; i < vect_tenor_Dates.size() ; ++i)
	{
		std::cout << "vect tenor dates i = " << i << "  " << vect_tenor_Dates[i] << std::endl ;
	} 

	std::vector<double> vect_delta_tau = simulationStructure->get_deltaT() ;
	for (size_t i = 0 ; i < vect_delta_tau.size() ; ++i)
	{
		std::cout << "vect delta tau i = " << i << "  " << vect_delta_tau[i] << std::endl ;
	} 
	VanillaSwaption_PTR vanillaSwaption_PTR_Test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;

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
	
	std::vector<LMM::Index> vect_fixedLegPaymentIndexSchedule, vect_floatLegPaymentIndexSchedule ;
	vect_fixedLegPaymentIndexSchedule = swap.get_fixedLegPaymentIndexSchedule();
	vect_floatLegPaymentIndexSchedule = swap.get_floatingLegPaymentIndexSchedule() ;
	std::cout << "indices des flux fixes" << std::endl ;
	for (size_t i= 0 ; i < vect_fixedLegPaymentIndexSchedule.size() ; ++i)
	{
		std::cout << vect_fixedLegPaymentIndexSchedule[i] << std::endl ;
	}		//retourne 2 4 6 8
	std::cout << "indices des flux flottants" << std::endl ;
	for (size_t i= 0 ; i < vect_floatLegPaymentIndexSchedule.size() ; ++i)
	{
		std::cout << vect_floatLegPaymentIndexSchedule[i] << std::endl ;
	}		//retourne 1 2 3 4 5 6 7 8

}


VanillaSwaption_PTR createSwaption()
{
	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 4; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 15) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	
	VanillaSwaption_PTR vanillaSwaption_PTR_Test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	vanillaSwaption_PTR_Test->show() ;

	return vanillaSwaption_PTR_Test ;
}

CheyetteDD_VanillaSwaptionApproxPricer_PTR createApproxPricer_PTR()
{
	CheyetteDD_Model_PTR	cheyetteDD_Model_Test_PTR(createCheyetteDD_Model()) ;
	VanillaSwaption_PTR		vanillaSwaption_Test_PTR(createSwaption()) ;
	CheyetteDD_VanillaSwaptionApproxPricer_PTR approxPricerTest_PTR(
							new CheyetteDD_VanillaSwaptionApproxPricer(	cheyetteDD_Model_Test_PTR, 
																		vanillaSwaption_Test_PTR)); 
	return approxPricerTest_PTR ;
}

void test_y_barre()
{
	std::vector<double> x, y, m_y ;
	x.push_back(0) ; x.push_back(1) ; x.push_back(2) ; 
	y.push_back(0.25) ; y.push_back(0.5) ;
	m_y.push_back(0) ; m_y.push_back(0) ;
	double k(1) ;

	int curveChoice = 1 ;  //0,1 ou 2
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	VanillaSwaption_PTR swaption_PTR_test(createSwaption()) ;

//cas m(t) = 0
	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;

	CheyetteDD_VanillaSwaptionApproxPricer approx = 
				CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, swaption_PTR_test); 

	double t = 1 ;
	double r0 = modele_test_PTR->get_courbeInput_PTR()->get_f_0_t(0) ;
	double res_integrale(exp(- 2) * pow(0.25 * r0, 2) * (exp(2) - 1) / 2) ;

	std::cout << "integrale_main   : " << res_integrale << std::endl ;
	std::cout << "integrale_classe : " << approx.get_buffer_y_bar_t(t) << std::endl ;

}

void test_ZC_swapRate_Num_Denom()
{
	std::vector<double> listeMatu, tauxZC ;
	listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100) ; 
	listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100) ; 
	listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100) ; 
	listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100) ;  
	listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100) ; 
	listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100) ; 
	listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100) ; 
	listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100) ;  
	listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100) ;
	listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100) ;
	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));

	std::vector<double> x, sigma_y, m_y ;
	x.push_back(0) ; x.push_back(1) ; x.push_back(2) ; 
	m_y.push_back(0.25) ; m_y.push_back(0.5) ;
	sigma_y.push_back(0.25) ; sigma_y.push_back(0.5) ;

	Piecewiseconst_RR_Function sigma = Piecewiseconst_RR_Function(x, sigma_y) ; 
	Piecewiseconst_RR_Function m = Piecewiseconst_RR_Function(x, m_y) ; 
	double k = 0.2 ;

	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR cheyetteDD_Model_PTR_Test(new CheyetteDD_Model( courbe_PTR_test, 
														CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m), shiftChoice)) ;
	cheyetteDD_Model_PTR_Test->show() ;


	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 15) );
	VanillaSwap swap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	
	VanillaSwaption_PTR vanillaSwaption_PTR_Test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	vanillaSwaption_PTR_Test->show() ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR approxPricerTest_PTR(
		new CheyetteDD_VanillaSwaptionApproxPricer(cheyetteDD_Model_PTR_Test, vanillaSwaption_PTR_Test)) ;
	
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
	double y_bar_t = approxPricerTest_PTR->get_buffer_y_bar_t(0.5) ;
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
	std::vector<double> listeMatu, tauxZC ;
	listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100) ; 
	listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100) ; 
	listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100) ; 
	listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100) ;  
	listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100) ; 
	listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100) ; 
	listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100) ; 
	listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100) ;  
	listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100) ;
	listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100) ;
	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));

	std::vector<double> x, sigma_y, m_y ;
	x.push_back(0) ; x.push_back(3) ; x.push_back(4) ; 
	m_y.push_back(0.2) ; m_y.push_back(0.2) ;
	sigma_y.push_back(0.25) ; sigma_y.push_back(0.25) ;

	Piecewiseconst_RR_Function sigma = Piecewiseconst_RR_Function(x, sigma_y) ; 
	Piecewiseconst_RR_Function m = Piecewiseconst_RR_Function(x, m_y) ; 
	double k = 0.2 ;

	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR cheyetteDD_Model_PTR_Test(new CheyetteDD_Model( courbe_PTR_test, 
													CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m), shiftChoice )) ;
	cheyetteDD_Model_PTR_Test->show() ;


	double strike          = 0.04;
	LMM::Index  indexStart = 4; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 5) );
	VanillaSwap swap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	
	VanillaSwaption_PTR vanillaSwaption_PTR_Test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	vanillaSwaption_PTR_Test->show() ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR approxPricerTest_PTR(
		new CheyetteDD_VanillaSwaptionApproxPricer(cheyetteDD_Model_PTR_Test, vanillaSwaption_PTR_Test)) ;

	std::cout << "approxPricerTest_PTR->swapRate(0, 0) : " << approxPricerTest_PTR->swapRate(0., 0., 0.) <<std::endl ;   //OK
	std::cout << "approxPricerTest_PTR->swapRate(0, 0) : " << approxPricerTest_PTR->swapRate(0.1, 0., 0.) <<std::endl ;  
}

void test_fonction_inverse()
{
	std::vector<double> x, y, m_y ;
	x.push_back(0) ; x.push_back(1) ; x.push_back(2) ; 
	y.push_back(0.25) ; y.push_back(0.5) ;
	m_y.push_back(0) ; m_y.push_back(0) ;
	double k(1) ;

	int curveChoice = 1 ;
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	VanillaSwaption_PTR swaption_PTR_test(createSwaption()) ;

	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;

	CheyetteDD_VanillaSwaptionApproxPricer approx = 
				CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, swaption_PTR_test); 

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
	std::vector<double> x, y, m_y ;
	x.push_back(0) ; x.push_back(1) ; x.push_back(2) ; 
	y.push_back(0.25) ; y.push_back(0.5) ;
	m_y.push_back(0) ; m_y.push_back(0) ;
	double k(1) ;

	CourbeInput_PTR courbe_PTR_test(createCourbeInput(1));
	VanillaSwaption_PTR swaption_PTR_test(createSwaption()) ;

	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;

	CheyetteDD_VanillaSwaptionApproxPricer approx = 
				CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, swaption_PTR_test); 

	std::cout << " " << std::endl ;
	std::cout << "------------  d sigma_r / dx  --------------" << std::endl ;
	std::cout << " " << std::endl ;

	std::cout << approx.get_CheyetteDD_Model()->sigma_r_t_1stDerivative(0., 0.5, 0.)  << std::endl ;
	std::cout << approx.get_CheyetteDD_Model()->sigma_r_t_1stDerivative(0.5, 0.5, 0.)  << std::endl ;
	std::cout << approx.get_CheyetteDD_Model()->sigma_r_t_1stDerivative(1., 0.5, 0.)  << std::endl ;

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
//swaption tres simplifiée
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(1));
	double strike          = 0.9/100;  //prendre ATM, proche de S0
	LMM::Index  indexStart = 4; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 5) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	swaption_PTR_test->show() ;

	std::vector<double> x, sigma_y, m_y ;
	x.push_back(0) ; x.push_back(3) ; x.push_back(4) ; 
//param
	double m_param = 0.0 ;
	double k = 0.2 ;
	double sigm = 0.35 ;

	m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	sigma_y.push_back(sigm) ; sigma_y.push_back(sigm) ;
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ; 


	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

	//valide pour m(t) = 0 uniquement
	//double t = 1 ;
	//std::cout	<< "y barre(t) should be " << sigm * sigm * 0.8/100 * 0.8/100 * (1 - exp(- 2 * k * t))/(2 * k) 
	//			<< " , is " << approx.get_buffer_y_bar_t(t)  
	//			<< std::endl ; 

	std::cout << "prixSwaption" << std::endl ;
	std::cout << approx.prixSwaptionApproxPiterbarg() << std::endl ;
	std::cout << "  " << std::endl ;

}
