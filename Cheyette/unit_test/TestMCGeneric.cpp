#include "TestMCGeneric.h"



//test swap maturité 5Y 6M contre fixe OK
//OK pour m(t) = 0, m(t) = 0.5 et m(t) = 1 (sigma = 20% et 50%)
void testgenericSwap()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	std::vector<double> x(2), sigma_y(1), m_y(1) ;
	x[0] = 0. ; x[1] = 10. ; 
	sigma_y[0] = 0.2 ;
	m_y[0] = 1. ;
	double k(0.2) ;

	std::vector<double> listeMatu, tauxZC ;
	listeMatu.push_back(0.) ;	tauxZC.push_back(0.8/100) ; 
	listeMatu.push_back(1.) ;	tauxZC.push_back(0.85/100) ; 
	listeMatu.push_back(2.) ;	tauxZC.push_back(0.9/100) ; 
	listeMatu.push_back(3.) ;	tauxZC.push_back(0.92/100) ;  
	listeMatu.push_back(4.) ;	tauxZC.push_back(0.95/100) ; 
	listeMatu.push_back(5.) ;	tauxZC.push_back(1.00/100) ; 
	listeMatu.push_back(10.) ;	tauxZC.push_back(1.5/100) ; 
	listeMatu.push_back(15.) ;	tauxZC.push_back(2.0/100) ;  
	listeMatu.push_back(20.) ;	tauxZC.push_back(2.5/100) ;
	listeMatu.push_back(25.) ;	tauxZC.push_back(2.3/100) ;
	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));

	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, sigma_y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct = CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;

	double fwdProbaT = 5. ;
	size_t indexMax = 10 ;
	std::vector<size_t>			simulationIndex(indexMax + 1), discretizationBetweenDates(indexMax + 1) ;	
	discretizationBetweenDates[0] = 0 ;
	simulationIndex[0] = 0 ;
	for (size_t i = 0 ; i <= indexMax ; ++i)
	{
		simulationIndex[i] = i ;
		discretizationBetweenDates[i] = 180 ; //180 dt env. = 1jour 
	}
	MC_Cheyette_PTR MC_Cheyette_PTR_Test( new MC_Cheyette(	modele_test_PTR,
															rnGenerator,
															Tenor::_6M,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates   ) ) ;

	MC_CheyetteDD_GenericSwapPricer_PTR mc(new MC_CheyetteDD_GenericSwapPricer(modele_test_PTR,
															rnGenerator,
															Tenor::_6M,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates )) ;
	
	double nominal(1.) ;
	bool ifFloored(false),	ifCapped(false);						
	double floorStrike(0.), capStrike(0.) ;
	double multiFactor(1.),addFactor(0.) ;
	double strike = 0.02 ;//0.009 ;
	Tenor floatingLegTenor = Tenor::_6M ;
	Tenor fixedLegTenor = Tenor::_12M ;
	LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, 5)) ;
	size_t indexStart(0), indexEnd(indexMax) ;

	//Instrument Factory construit un swap payeur (floating - fixed)
	GeneticSwap_CONSTPTR genericSwapTest(InstrumentFactory::createVanillaSwap(	strike, 
																				indexStart, indexEnd, 
																				floatingLegTenor, fixedLegTenor,
																				swapStructure, nominal));

//comparaison avec Instrument Factory OK :

//	size_t valuationDateIndex(0) ;
////float leg	
//	std::vector<Coupon_CONSTPTR>  legFloat ;
//	double deltaFloat = 0.5 ;
//
//	for (size_t i = 1 ; i <= indexMax ; ++i)
//	{
//		size_t fixingTime	= i - 1 ;
//		size_t paymentIndex = i ;
//		Rate_CONSTPTR rate (new LiborRate(fixingTime, Tenor::_6M)) ;
//		Coupon_CONSTPTR couponFloat(new CappedFlooredCoupon(paymentIndex, nominal, deltaFloat, ifFloored, floorStrike,
//													ifCapped, capStrike, rate, multiFactor, addFactor, valuationDateIndex) );	
//		legFloat.push_back(couponFloat) ;
//	}
//
////fixed leg
//	std::vector<Coupon_CONSTPTR>  legFixed ;
//	double deltaFixed = 1. ;
//	
//	for (size_t i = 2 ; i <= indexMax ; i+=2)		//à modifier selon les tenors fixed et float
//	{
//		size_t paymentIndex = i ;
//		Rate_CONSTPTR rate (new ConstRate(strike)) ;
//		Coupon_CONSTPTR couponFixed(new CappedFlooredCoupon(paymentIndex, nominal, deltaFixed, ifFloored, floorStrike,
//													ifCapped, capStrike, rate, multiFactor, addFactor, valuationDateIndex) );
//		legFixed.push_back(couponFixed) ;
//	}
//
//	//swap payeur
//	CouponLeg_CONSTPTR couponLeg1(new CouponLeg(legFloat));	
//	CouponLeg_CONSTPTR couponLeg2(new CouponLeg(legFixed));	
//
//	GeneticSwap_CONSTPTR 	genericSwapTest(new GeneticSwap(couponLeg1, couponLeg2)) ;

	//std::cout << "prix swap avec courbe input : " <<	strike * exp(- 2 * 0.9/100)
	//												-	(exp(- 0.85/100) - exp(- 2 * 0.9/100)) << std::endl ;


	double fixedL = 0. ;
	fixedL += 1. * exp(- 1. * 0.85/100) ;
	fixedL += 1. * exp(- 2. * 0.9/100) ;
	fixedL += 1. * exp(- 3. * 0.92/100) ;
	fixedL += 1. * exp(- 4. * 0.95/100) ;
	fixedL += 1. * exp(- 5. * 1.0/100) ;
	fixedL = strike * fixedL ;
	double floatL = 1 - exp(- 5 * 1.00/100) ;
	std::cout << "prix swap avec courbe input : " << floatL - fixedL << std::endl ;

	size_t nbSimulation = 10000 ;
	mc->swapNPV(genericSwapTest, nbSimulation, floatingLegTenor, fixedLegTenor) ;
	nbSimulation = 40000 ;
	mc->swapNPV(genericSwapTest, nbSimulation, floatingLegTenor, fixedLegTenor) ;
	nbSimulation = 100000 ;
	mc->swapNPV(genericSwapTest, nbSimulation, floatingLegTenor, fixedLegTenor) ;
}

//OK
void testgenericSwaption()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	std::vector<double> x(2), sigma_y(1), m_y(1) ;
	x[0] = 0. ; x[1] = 10. ; 
	sigma_y[0] = 0.2 ;
	m_y[0] = 1. ;
	double k(0.2) ;
	
	std::vector<double> listeMatu, tauxZC ;
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

	listeMatu.push_back(0.) ;	tauxZC.push_back(0.8/100) ; 
	listeMatu.push_back(1.) ;	tauxZC.push_back(0.85/100) ; 
	listeMatu.push_back(2.) ;	tauxZC.push_back(0.9/100) ; 
	listeMatu.push_back(3.) ;	tauxZC.push_back(0.92/100) ;  
	listeMatu.push_back(4.) ;	tauxZC.push_back(0.95/100) ; 
	listeMatu.push_back(5.) ;	tauxZC.push_back(1.00/100) ; 
	listeMatu.push_back(10.) ;	tauxZC.push_back(1.5/100) ; 
	listeMatu.push_back(15.) ;	tauxZC.push_back(2.0/100) ;  
	listeMatu.push_back(20.) ;	tauxZC.push_back(2.5/100) ;
	listeMatu.push_back(25.) ;	tauxZC.push_back(2.3/100) ;
	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));

	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

	double fwdProbaT = 5. ;
	size_t indexMax = 10 ;
	std::vector<size_t>			simulationIndex(indexMax + 1), discretizationBetweenDates(indexMax + 1) ;	
	discretizationBetweenDates[0] = 0 ;
	simulationIndex[0] = 0 ;
	for (size_t i = 0 ; i <= indexMax ; ++i)
	{
		simulationIndex[i] = i ;
		discretizationBetweenDates[i] = 20 ; //180 dt env. = 1jour 
	}

	MC_CheyetteDD_GenericSwaptionPricer_PTR mc(new MC_CheyetteDD_GenericSwaptionPricer(modele_test_PTR,
															rnGenerator,
															Tenor::_6M,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates   ) ) ;

	double nominal(1.) ;
	bool ifFloored(false),	ifCapped(false);						
	double floorStrike(0.), capStrike(0.) ;
	double multiFactor(1.),addFactor(0.) ;
	double strike = -10 ;//0.009 ;
	Tenor floatingLegTenor = Tenor::_6M ;
	Tenor fixedLegTenor = Tenor::_12M ;
	LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, 5)) ;
	size_t indexStart(0), indexEnd(indexMax) ;

	//Instrument Factory construit un swap payeur (floating - fixed)
	GeneticSwap_CONSTPTR genericSwapTest(InstrumentFactory::createVanillaSwap(	strike, 
																				indexStart, indexEnd, 
																				floatingLegTenor, fixedLegTenor,
																				swapStructure, nominal));

	size_t maturity = 4 ; //index
	GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;

	size_t nbMC = 3 ;
	std::vector<size_t> nbSimus(nbMC) ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	nbSimus[0] = 10000 ; nbSimus[1] = 30000 ; nbSimus[2] = 50000 ;

	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
	}
	
	//std::cout << "prix swap avec courbe input : " <<		(exp(- 2 * 0.9/100) - exp(- 3 * 0.92/100)) 
	//												-  strike * exp(- 3 * 0.92/100) << std::endl ;

	double fixedL = 0. ;
	fixedL += 1. * exp(- 1. * 0.85/100) ;
	fixedL += 1. * exp(- 2. * 0.9/100) ;
	fixedL += 1. * exp(- 3. * 0.92/100) ;
	fixedL += 1. * exp(- 4. * 0.95/100) ;
	fixedL += 1. * exp(- 5. * 1.0/100) ;
	fixedL = strike * fixedL ;
	double floatL = 1 - exp(- 5 * 1.00/100) ;
	std::cout << "prix swap avec courbe input : " << floatL - fixedL << std::endl ;

	mc->print(genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
}


void testVanillaSwap()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	std::vector<double> x(2), sigma_y(1), m_y(1) ;
	x[0] = 0. ; x[1] = 10. ; 
	sigma_y[0] = 0.2 ;
	m_y[0] = 1. ;
	double k(0.2) ;

	std::vector<double> listeMatu, tauxZC ;
	listeMatu.push_back(0.) ;	tauxZC.push_back(0.8/100) ; 
	listeMatu.push_back(1.) ;	tauxZC.push_back(0.85/100) ; 
	listeMatu.push_back(2.) ;	tauxZC.push_back(0.9/100) ; 
	listeMatu.push_back(3.) ;	tauxZC.push_back(0.92/100) ;  
	listeMatu.push_back(4.) ;	tauxZC.push_back(0.95/100) ; 
	listeMatu.push_back(5.) ;	tauxZC.push_back(1.00/100) ; 
	listeMatu.push_back(10.) ;	tauxZC.push_back(1.5/100) ; 
	listeMatu.push_back(15.) ;	tauxZC.push_back(2.0/100) ;  
	listeMatu.push_back(20.) ;	tauxZC.push_back(2.5/100) ;
	listeMatu.push_back(25.) ;	tauxZC.push_back(2.3/100) ;
	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));

	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, sigma_y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct = CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;

	double fwdProbaT = 5. ;
	size_t indexMax = 10 ;
	std::vector<size_t>			simulationIndex(indexMax + 1), discretizationBetweenDates(indexMax + 1) ;	
	discretizationBetweenDates[0] = 0 ;
	simulationIndex[0] = 0 ;
	for (size_t i = 0 ; i <= indexMax ; ++i)
	{
		simulationIndex[i] = i ;
		discretizationBetweenDates[i] = 180 ; //180 dt env. = 1jour 
	}
	MC_Cheyette_PTR MC_Cheyette_PTR_Test( new MC_Cheyette(	modele_test_PTR,
															rnGenerator,
															Tenor::_6M,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates   ) ) ;

	MC_CheyetteDD_VanillaSwapPricer_PTR mc(new MC_CheyetteDD_VanillaSwapPricer(modele_test_PTR,
															rnGenerator,
															Tenor::_6M,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates )) ;
	
	double strike = 0.02 ;//0.009 ;
	Tenor floatingLegTenor = Tenor::_6M ;
	Tenor fixedLegTenor = Tenor::_12M ;
	LMMTenorStructure_PTR swapStructure(new LMMTenorStructure(floatingLegTenor, 5)) ;
	size_t indexStart(0), indexEnd(indexMax) ;

		//VanillaSwap(const double& strike,
		//			LMM::Index  indexStart, 
		//			LMM::Index  indexEnd, 
		//			const Tenor& floatingLegTenorType,		
		//			const Tenor& fixedLegTenorType,
		//			LMMTenorStructure_PTR lmmTenorStructure); 


	VanillaSwap_PTR vanillaSwapTest(new VanillaSwap(strike, indexStart, indexEnd, 
													floatingLegTenor, fixedLegTenor,
													swapStructure));

//	size_t valuationDateIndex(0) ;
////float leg	
//	std::vector<Coupon_CONSTPTR>  legFloat ;
//	double deltaFloat = 0.5 ;
//
//	for (size_t i = 1 ; i <= indexMax ; ++i)
//	{
//		size_t fixingTime	= i - 1 ;
//		size_t paymentIndex = i ;
//		Rate_CONSTPTR rate (new LiborRate(fixingTime, Tenor::_6M)) ;
//		Coupon_CONSTPTR couponFloat(new CappedFlooredCoupon(paymentIndex, nominal, deltaFloat, ifFloored, floorStrike,
//													ifCapped, capStrike, rate, multiFactor, addFactor, valuationDateIndex) );	
//		legFloat.push_back(couponFloat) ;
//	}
//
////fixed leg
//	std::vector<Coupon_CONSTPTR>  legFixed ;
//	double deltaFixed = 1. ;
//	
//	for (size_t i = 2 ; i <= indexMax ; i+=2)		//à modifier selon les tenors fixed et float
//	{
//		size_t paymentIndex = i ;
//		Rate_CONSTPTR rate (new ConstRate(strike)) ;
//		Coupon_CONSTPTR couponFixed(new CappedFlooredCoupon(paymentIndex, nominal, deltaFixed, ifFloored, floorStrike,
//													ifCapped, capStrike, rate, multiFactor, addFactor, valuationDateIndex) );
//		legFixed.push_back(couponFixed) ;
//	}
//
//	//swap payeur
//	CouponLeg_CONSTPTR couponLeg1(new CouponLeg(legFloat));	
//	CouponLeg_CONSTPTR couponLeg2(new CouponLeg(legFixed));	
//
//	GeneticSwap_CONSTPTR 	genericSwapTest(new GeneticSwap(couponLeg1, couponLeg2)) ;

	//std::cout << "prix swap avec courbe input : " <<	strike * exp(- 2 * 0.9/100)
	//												-	(exp(- 0.85/100) - exp(- 2 * 0.9/100)) << std::endl ;


	double fixedL = 0. ;
	fixedL += 1. * exp(- 1. * 0.85/100) ;
	fixedL += 1. * exp(- 2. * 0.9/100) ;
	fixedL += 1. * exp(- 3. * 0.92/100) ;
	fixedL += 1. * exp(- 4. * 0.95/100) ;
	fixedL += 1. * exp(- 5. * 1.0/100) ;
	fixedL = strike * fixedL ;
	double floatL = 1 - exp(- 5 * 1.00/100) ;
	std::cout << "prix swap avec courbe input : " << floatL - fixedL << std::endl ;

	size_t nbSimulation = 10000 ;
	mc->swapNPV(vanillaSwapTest, nbSimulation, floatingLegTenor, fixedLegTenor) ;
	nbSimulation = 40000 ;
	mc->swapNPV(vanillaSwapTest, nbSimulation, floatingLegTenor, fixedLegTenor) ;
	nbSimulation = 100000 ;
	mc->swapNPV(vanillaSwapTest, nbSimulation, floatingLegTenor, fixedLegTenor) ;
}

void testVanillaSwaption()
{

}

