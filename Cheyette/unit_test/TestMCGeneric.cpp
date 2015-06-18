#include "TestMCGeneric.h"




void testgenericSwap()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	std::vector<double> x, y ;
	x.push_back(0) ; x.push_back(1) ; x.push_back(2) ; 
	y.push_back(0.2) ; y.push_back(0.25) ;
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

	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct = CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m) ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct)) ;

	double fwdProbaT = 2. ;
	std::vector<size_t>			simulationIndex ;	
	simulationIndex.push_back(0) ; simulationIndex.push_back(2) ; simulationIndex.push_back(3) ; simulationIndex.push_back(4) ;

	std::vector<size_t>			discretizationBetweenDates ;
	discretizationBetweenDates.push_back(0) ; 
	discretizationBetweenDates.push_back(100) ; 
	discretizationBetweenDates.push_back(50) ;
	discretizationBetweenDates.push_back(50) ;

	MC_Cheyette_PTR MC_Cheyette_PTR_Test( new MC_Cheyette(	modele_test_PTR,
															rnGenerator,
															Tenor::_6M,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates   ) ) ;

	MC_CheyetteDD_GenericSwapPricer_PTR mc(new MC_CheyetteDD_GenericSwapPricer(MC_Cheyette_PTR_Test)) ;
	
	double nominal(1) ;
	bool ifFloored(false),	ifCapped(false);						
	double floorStrike(0), capStrike(0) ;
	double multiFactor(1),addFactor(0) ;
	size_t valuationDateIndex(0) ;

	size_t paymentIndex(4) ;
	double deltaFixed = 1 ;
	double strike = 0.02 ;//0.009 ;
	Rate_CONSTPTR rate (new ConstRate(strike)) ;
	Coupon_CONSTPTR couponFixed(new CappedFlooredCoupon(paymentIndex, nominal, deltaFixed, ifFloored, floorStrike,
												ifCapped, capStrike, rate, multiFactor, addFactor, valuationDateIndex) );

	paymentIndex = 3 ;
	double deltaFloat = 0.5 ;
	size_t fixingTime(2) ;
	Rate_CONSTPTR rate1 (new LiborRate(fixingTime, Tenor::_6M)) ;
	Coupon_CONSTPTR couponFloat1(new CappedFlooredCoupon(paymentIndex, nominal, deltaFloat, ifFloored, floorStrike,
												ifCapped, capStrike, rate1, multiFactor, addFactor, valuationDateIndex) );

	paymentIndex = 4 ;
	fixingTime = 3  ;
	Rate_CONSTPTR rate2 (new LiborRate(fixingTime, Tenor::_6M)) ;
	Coupon_CONSTPTR couponFloat2(new CappedFlooredCoupon(paymentIndex, nominal, deltaFloat, ifFloored, floorStrike,
												ifCapped, capStrike, rate2, multiFactor, addFactor, valuationDateIndex) );
	
	std::vector<Coupon_CONSTPTR>  leg1 ;
	leg1.push_back(couponFixed) ;
	std::vector<Coupon_CONSTPTR>  leg2 ;
	leg2.push_back(couponFloat1) ; leg2.push_back(couponFloat2) ;
	CouponLeg_CONSTPTR couponLeg1(new CouponLeg(leg1));
	CouponLeg_CONSTPTR couponLeg2(new CouponLeg(leg2));

	GeneticSwap_CONSTPTR 	geneticSwapTest(new GeneticSwap(couponLeg1, couponLeg2)) ;

	std::cout << "prix swap avec courbe input : " <<	strike * exp(- 2 * 0.9/100)
													-	(exp(- 0.85/100) - exp(- 2 * 0.9/100)) << std::endl ;
	//size_t nbSimulation = 5000 ;
	//std::cout << "prix MC swap :                " << mc->swapNPV(geneticSwapTest, nbSimulation) << std::endl ;
	//nbSimulation = 10000 ;
	//std::cout << "prix MC swap :                " << mc->swapNPV(geneticSwapTest, nbSimulation) << std::endl ;
	//nbSimulation = 20000 ;
	//std::cout << "prix MC swap :                " << mc->swapNPV(geneticSwapTest, nbSimulation) << std::endl ;

//	size_t nbSimulation = 5000 ;
	//mc->swapNPV(geneticSwapTest, nbSimulation) ;
	//nbSimulation = 10000 ;
	//mc->swapNPV(geneticSwapTest, nbSimulation) ;
	//nbSimulation = 20000 ;
	//mc->swapNPV(geneticSwapTest, nbSimulation) ;
	size_t	nbSimulation = 100000 ;
	mc->swapNPV(geneticSwapTest, nbSimulation) ;
	nbSimulation = 200000 ;
	mc->swapNPV(geneticSwapTest, nbSimulation) ;
}


void testgenericSwaption()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	std::vector<double> x, m_y, sigma_y ;
	x.push_back(0) ; x.push_back(3) ; x.push_back(4) ; 
//param
	double m_param = 0.0 ;
	double k = 0.2 ;
	double sigm = 0.35 ;

	m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	sigma_y.push_back(sigm) ; sigma_y.push_back(sigm) ;
	
	std::vector<double> listeMatu, tauxZC ;


	double translation = 0.0 ;
//courbe des taux plate
	std::cout << "courbe des taux plate à 1 %" << std::endl ;
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

	//listeMatu.push_back(0.) ;	tauxZC.push_back(0.8/100) ; 
	//listeMatu.push_back(1.) ;	tauxZC.push_back(0.85/100) ; 
	//listeMatu.push_back(2.) ;	tauxZC.push_back(0.9/100) ; 
	//listeMatu.push_back(3.) ;	tauxZC.push_back(0.92/100) ;  
	//listeMatu.push_back(4.) ;	tauxZC.push_back(0.95/100) ; 
	//listeMatu.push_back(5.) ;	tauxZC.push_back(1.00/100) ; 
	//listeMatu.push_back(10.) ;	tauxZC.push_back(1.5/100) ; 
	//listeMatu.push_back(15.) ;	tauxZC.push_back(2.0/100) ;  
	//listeMatu.push_back(20.) ;	tauxZC.push_back(2.5/100) ;
	//listeMatu.push_back(25.) ;	tauxZC.push_back(2.3/100) ;
	//CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));

	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct)) ;
	modele_test_PTR->show() ;

	double fwdProbaT = 3. ;
	std::vector<size_t>			simulationIndex ;	
	simulationIndex.push_back(0) ; 
	simulationIndex.push_back(4) ; 
	simulationIndex.push_back(5) ; 
	simulationIndex.push_back(6) ;

	std::vector<size_t>			discretizationBetweenDates ;
	discretizationBetweenDates.push_back(0) ; 
	discretizationBetweenDates.push_back(200) ; 
	discretizationBetweenDates.push_back(50) ;
	discretizationBetweenDates.push_back(50) ;

	MC_Cheyette_PTR MC_Cheyette_PTR_Test( new MC_Cheyette(	modele_test_PTR,
															rnGenerator,
															Tenor::_6M,
															fwdProbaT,
															simulationIndex,		
															discretizationBetweenDates   ) ) ;

	MC_CheyetteDD_GenericSwaptionPricer_PTR mc(new MC_CheyetteDD_GenericSwaptionPricer(MC_Cheyette_PTR_Test)) ;
	
	double nominal(1) ;
	bool ifFloored(false),	ifCapped(false);						
	double floorStrike(0), capStrike(0) ;
	double multiFactor(1),addFactor(0) ;
	size_t valuationDateIndex(0) ;


	size_t paymentIndex(5) ;
	double deltaFloat = 0.5 ;
	size_t fixingTime(4) ;
	Rate_CONSTPTR rate1 (new LiborRate(fixingTime, Tenor::_6M)) ;
	Coupon_CONSTPTR couponFloat1(new CappedFlooredCoupon(paymentIndex, nominal, deltaFloat, ifFloored, floorStrike,
												ifCapped, capStrike, rate1, multiFactor, addFactor, valuationDateIndex) );

	paymentIndex = 6 ;
	fixingTime = 5  ;
	Rate_CONSTPTR rate2 (new LiborRate(fixingTime, Tenor::_6M)) ;
	Coupon_CONSTPTR couponFloat2(new CappedFlooredCoupon(paymentIndex, nominal, deltaFloat, ifFloored, floorStrike,
												ifCapped, capStrike, rate2, multiFactor, addFactor, valuationDateIndex) );

	paymentIndex = 6 ;
	double deltaFixed = 1 ;
	double strike = 0.9/100 ;
	Rate_CONSTPTR rate (new ConstRate(strike)) ;
	Coupon_CONSTPTR couponFixed(new CappedFlooredCoupon(paymentIndex, nominal, deltaFixed, ifFloored, floorStrike,
												ifCapped, capStrike, rate, multiFactor, addFactor, valuationDateIndex) );

	std::vector<Coupon_CONSTPTR>  leg1 ;
	leg1.push_back(couponFloat1) ; 
	leg1.push_back(couponFloat2) ;

	std::vector<Coupon_CONSTPTR>  leg2 ;
	leg2.push_back(couponFixed) ;

	CouponLeg_CONSTPTR couponLeg1(new CouponLeg(leg1));
	CouponLeg_CONSTPTR couponLeg2(new CouponLeg(leg2));

	GeneticSwap_CONSTPTR 	genericSwapTest(new GeneticSwap(couponLeg1, couponLeg2)) ; 

	size_t maturity = 4 ; //index
	GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;

	size_t nbMC = 3 ;
	std::vector<size_t> nbSimus(nbMC) ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	nbSimus[0] = 10 ; nbSimus[1] = 15 ; nbSimus[2] = 30 ;

	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
	}
	
	//std::cout << "prix swap avec courbe input : " <<		(exp(- 2 * 0.9/100) - exp(- 3 * 0.92/100)) 
	//												-  strike * exp(- 3 * 0.92/100) << std::endl ;

	mc->print(genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
}