#include "Test_comparaison_MC_approx.h"

//comparaison prix approx et MC
//Modele Cheyette DD ou Quad, approx lineaire ou quadratique
void test_approx()
{	
	size_t xmax = 4 ;
	int numCourbe = 1 ;
	double k = 0.02 ;
	
	double sigmaValue = 0.2 ;
	double mValue = 0. ;
	CheyetteDD_Model_PTR		cheyetteModel_PTR = creeCheyetteDD_Modele_PTR(xmax, numCourbe, k, sigmaValue, mValue) ;

	//double aValue = 0.2 ;
	//double bValue = 0.2 ;
	//double cValue = 0.2 ;
	//CheyetteQuad_Model_PTR	cheyetteModel_PTR = creeCheyetteQuad_Modele_PTR(xmax, numCourbe, k, aValue, bValue, cValue) ;

//swaption
	VanillaSwaption_PTR pVanillaSwaption = createSwaptionTest() ;    //strike = 0.04;

//approx
	Cheyette_SwaptionPricer_LinearApprox_PTR approxPricerTest_PTR(
							new Cheyette_SwaptionPricer_LinearApprox(cheyetteModel_PTR)  ) ; 

	//Cheyette_SwaptionPricer_QuadApprox_PTR approxPricerTest_PTR(
	//						new Cheyette_SwaptionPricer_QuadApprox(cheyetteModel_PTR)  ) ; 

	double prixApprox = approxPricerTest_PTR->price(pVanillaSwaption) ;
	std::cout << "Prix swaption approx : " << prixApprox << std::endl ;

//MC
	LMMTenorStructure_PTR	pTenorStructure(new LMMTenorStructure(Tenor::_6M, static_cast<int>(xmax))) ;
	size_t					fwdProbaT = static_cast<size_t>(xmax) ;
	size_t					discretizationBetweenDates = 200 ;

	MC_Cheyette_VanillaSwaptionPricer_PTR pSwaptionPricer = creeMC_SwaptionPricer_PTR(	cheyetteModel_PTR, 
																				pTenorStructure,
																				fwdProbaT,
																				discretizationBetweenDates) ;

	size_t nbSimu = 300000 ;
	double MCswaptionPrice = pSwaptionPricer->price(pVanillaSwaption, nbSimu)[0] ;
	std::cout << "Prix swaption MC : " << MCswaptionPrice << std::endl ;




//	double strike =  - 1000. ; //temporaire avant de mettre strike ATM
//	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
//	
//	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//	CheyetteDD_VanillaSwaptionApproxPricer approx(pCheyetteDD_Model, swaption_PTR_test);
//
////calcul du strike ATM pour le swap
//	double strikeATM = approx.swapRate0() ;
//
//	o << "strike ATM pour swaption " << a << "Y" << b << "Y : " << strikeATM << std::endl ;
//	std::cout << "strike ATM pour swaption " << a << "Y" << b << "Y : " << strikeATM << std::endl ;
//	std::cout << std::endl ;
//
//	swap.set_strike(strikeATM) ;
//	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//	swaption_PTR_test_ATM->show() ;
//	CheyetteDD_VanillaSwaptionApproxPricer approxATM(pCheyetteDD_Model, swaption_PTR_test_ATM);
//
//	std::cout << "prixSwaption" << std::endl ;
//	double approxPrice = approxATM.prixSwaptionApproxPiterbarg() ;
//	double b_barre = approxATM.get_buffer_b_barre_() ;
//	std::cout << "approxPrice : " << approxPrice << std::endl ;
//	std::cout << "b barre :     " << b_barre << std::endl ;
//	std::cout << "  " << std::endl ;
//
//
//	size_t nbMC = nbSimus.size() ;
//	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
//	for (size_t i = 0 ; i < nbMC ; ++i)
//	{
//		std::vector<double> resMC = mc->price(swaption_PTR_test_ATM, nbSimus[i]) ;
//		vectPrixMC[i] = resMC[0] ;
//		vectICinf[i]  = resMC[1] ;
//		vectICsup[i]  = resMC[2] ;
//	}
//	
//	double annuityA0 = approxATM.swapRateDenominator(0., 0., 0.) ;
//	double swapRateS0 = approxATM.swapRate0() ;
//
//	double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
//																		swapRateS0, strikeATM, a) ;
//
//	std::cout << "strike :    " << strikeATM << std::endl ;
//	std::cout << "expiry :    " << a << std::endl ;
//	std::cout << "annuity :   " << annuityA0 << std::endl ;
//	std::cout << "swapRate :  " << swapRateS0 << std::endl ;
//	std::cout << "vol Black : " << volBlack << std::endl ;
//
//	mc->printMC_vs_approx(o, approxPrice, b_barre, annuityA0, swapRateS0, volBlack,
//							a, b, swaption_PTR_test_ATM, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
}
