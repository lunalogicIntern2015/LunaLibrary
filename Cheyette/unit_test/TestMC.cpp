#include "TestMC.h"

#include <fstream>
#include <ostream>


MC_Cheyette_VanillaSwapPricer_PTR creeMC_SwapPricer_PTR(CheyetteModel_PTR		cheyetteModel_PTR, 
														LMMTenorStructure_PTR	pTenorStructure,
														size_t					fwdProbaT,
														size_t					discretizationBetweenDates)
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	MC_Cheyette_VanillaSwapPricer_PTR mcCheyette(new MC_Cheyette_VanillaSwapPricer(cheyetteModel_PTR,
																				rnGenerator, pTenorStructure, fwdProbaT,
																				discretizationBetweenDates)) ;
	return mcCheyette ;
}

MC_Cheyette_VanillaSwaptionPricer_PTR creeMC_SwaptionPricer_PTR(CheyetteModel_PTR		cheyetteModel_PTR, 
															LMMTenorStructure_PTR	pTenorStructure,
															size_t					fwdProbaT,
															size_t					discretizationBetweenDates)
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	MC_Cheyette_VanillaSwaptionPricer_PTR mcCheyette(new MC_Cheyette_VanillaSwaptionPricer(cheyetteModel_PTR,
																				rnGenerator, pTenorStructure, fwdProbaT,
																				discretizationBetweenDates)) ;
	return mcCheyette ;
}

//proba forward QT
void testMCSwap_SwaptionPricer()
{
	size_t xmax = 3 ;
	int numCourbe = 1 ;
	double k = 0.02 ;
	double sigmaValue = 0.2 ;
	double mValue = 0. ;

	LMMTenorStructure_PTR	pTenorStructure(new LMMTenorStructure(Tenor::_6M, static_cast<int>(xmax))) ;
	size_t					fwdProbaT = static_cast<size_t>(xmax) ;
	size_t					discretizationBetweenDates = 200 ;

	CheyetteModel_PTR		cheyetteModel_PTR = creeCheyetteDD_Modele_PTR(xmax, numCourbe, k, sigmaValue, mValue) ;
	MC_Cheyette_VanillaSwapPricer_PTR pSwapPricer = creeMC_SwapPricer_PTR(	cheyetteModel_PTR, 
																				pTenorStructure,
																				fwdProbaT,
																				discretizationBetweenDates) ;
	VanillaSwap_PTR pVanillaSwap = createSwapTest() ;
	VanillaSwaption_PTR pVanillaSwaption = createSwaptionTest() ;    //strike = 0.04;

	//float - fixed
	CourbeInput_PTR pCourbe = cheyetteModel_PTR->get_courbeInput_PTR() ;

	size_t indexStart = pVanillaSwap->get_indexStart() ;
	size_t indexEnd = pVanillaSwap->get_indexEnd() ;
	double tauxStart	= pCourbe->get_tauxZC0(indexStart) ;
	double tauxEnd		= pCourbe->get_tauxZC0(indexEnd / 2.) ;  //...

	double ZC_T0 = exp( - static_cast<double>(indexStart / 2.) * tauxStart ) ; 
	double ZC_TN = exp( - static_cast<double>((indexStart + indexEnd)/2.)  * tauxEnd ) ;
	double floatLeg = ZC_T0 - ZC_TN ;
	
	std::vector<size_t> fixedLegIndex = pVanillaSwap->get_fixedLegPaymentIndexSchedule() ;
	double fixedLeg = 0. ;
	for (size_t i = 0 ; i < fixedLegIndex.size() ; ++i)
	{
		double date	= fixedLegIndex[i] / 2. ;
		fixedLeg += exp( - date * pCourbe->get_tauxZC0(date) ) ;
	}

	std::cout << "Prix swap courbe : " << floatLeg - pVanillaSwap->get_strike() * fixedLeg << std::endl ;

//prix swap MC
	size_t nbSimu = 10000 ;
	double MCswapPrice = pSwapPricer->swapNPV(pVanillaSwap, nbSimu) ;
	std::cout << "Prix swap MC : " << MCswapPrice << std::endl ;

//swaption
	MC_Cheyette_VanillaSwaptionPricer_PTR pSwaptionPricer = creeMC_SwaptionPricer_PTR(	cheyetteModel_PTR, 
																				pTenorStructure,
																				fwdProbaT,
																				discretizationBetweenDates) ;

	double MCswaptionPrice = pSwaptionPricer->price(pVanillaSwaption, nbSimu)[0] ;
	std::cout << "Prix swaption MC : " << MCswaptionPrice << std::endl ;

	pVanillaSwap->set_strike(-10.) ;
	pVanillaSwaption->getUnderlyingSwap_RefNonConst().set_strike(-10.) ;

	MCswapPrice = pSwapPricer->swapNPV(pVanillaSwap, nbSimu) ;
	std::cout << "Prix swap MC : " << MCswapPrice << std::endl ;

	MCswaptionPrice = pSwaptionPricer->price(pVanillaSwaption, nbSimu)[0] ;
	std::cout << "Prix swaption MC : " << MCswaptionPrice << std::endl ;

}

//void testMCSwapPricer_annuity(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor)
//{
//	double strike          = 0.02 ;//- 1000. ;
//
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = creeModele(a, b, floatTenor, fixedTenor, strike) ;
//	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
//	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
//	VanillaSwap_PTR pSwap(new VanillaSwap(pSwaption->getUnderlyingSwap())) ;
//
//	double S0 = pApproxPricer->swapRate0() ;	
//	double r0 = pApproxPricer->get_CheyetteDD_Model()->get_courbeInput_PTR()->get_f_0_t(0) ;
//
//	double sigma = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().sigma_(0.5) ;//cste
//	double k = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().k_ ;
//	
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	
//	std::vector<double> res(3) ;
//	double somme_xi_swap	  = 0. ;
//	double somme_xi2_swap	  = 0. ;
//	double somme_xi_swaption  = 0. ;
//	double somme_xi2_swaption = 0. ;
//	
//	size_t nbPas = 200 ;
//	double dt = double(a)/double(nbPas) ;
//	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
//	{
//		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
//		std::vector<double>   gaussian_tmp(nbPas);  
//		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
//		double t = 0. ;
//		double S_t = S0 ;
//
//		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
//		{
//			double y_t, inv ;
//			if (t == 0.)
//			{
//				y_t = 0. ;
//				inv = 0. ;		
//			}
//			else 
//			{
//				y_t = sigma * sigma * r0 * r0 / (2 * k) * (1 - exp(-2 * k *t)) ;
//					pApproxPricer->get_buffer_y_bar_t(t) ;   // ...
//				inv = pApproxPricer->inverse(t, S_t) ;
//			}
//
//			double derivee = pApproxPricer->swapRate_1stDerivative(t, inv, y_t) ;
//			double sigma = pApproxPricer->get_CheyetteDD_Model()->sigma_r(t, inv, y_t) ;
//
//			double S_t_plus_dt = S_t + derivee * sigma * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
//			S_t = S_t_plus_dt ;
//			t += dt ;
//		}
//		double resSwap = S_t - strike ;
//		double resSwaption = std::max(S_t - strike, 0.) ; //  (S(T0) - K)+ 
//		
//		somme_xi_swap += resSwap ;
//		somme_xi2_swap += resSwap * resSwap ;
//
//		somme_xi_swaption += resSwaption ;
//		somme_xi2_swaption += resSwaption * resSwaption ;
//	}
//
//	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//swap
//	double mean_x_swap	= A0 * somme_xi_swap / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
//	double mean_x2_swap	= A0 * A0 * somme_xi2_swap / nbSimus; 
//	double variance_swap = mean_x2_swap - mean_x_swap * mean_x_swap ;
//
//	std::cout <<"variance : " << variance_swap << std::endl ;
//
//	double IC_left_swap	= mean_x_swap - 2.57*std::sqrt(variance_swap / nbSimus);
//	double IC_right_swap = mean_x_swap + 2.57*std::sqrt(variance_swap / nbSimus);
//
//	std::cout   << "SWAP " << std::endl;
//	std::cout   << "MC swap : " << mean_x_swap << std::endl;
//	std::cout	<< "nb simulations : " << nbSimus << std::endl;
//	std::cout   << "99% confidence interval  [" << IC_left_swap << " , " << IC_right_swap	<< "]" << std::endl;
//	
//	CourbeInput_PTR pCourbe = pApproxPricer->get_CheyetteDD_Model()->get_courbeInput_PTR() ;
//
//	double tauxStart	= pCourbe->get_tauxZC0(a) ;
//	double tauxEnd		= pCourbe->get_tauxZC0(a + b) ;
//	double ZC_T0 = exp( - double(a) * tauxStart ) ; //exp( - indexStart / 2. * tauxStart ) ;
//	double ZC_TN = exp( - double((a+b)) * tauxEnd ) ;
//	double floatLeg = ZC_T0 - ZC_TN ;
//	
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	std::vector<size_t> fixedLegIndex = pSwap->get_fixedLegPaymentIndexSchedule() ;
//	double fixedLeg = 0. ;
//	for (size_t i = 0 ; i < fixedLegIndex.size() ; ++i)
//	{
//		double date	= fixedLegIndex[i] * tenor ;
//		fixedLeg += exp( - date * pCourbe->get_tauxZC0(date) ) ;
//	}
//
//	std::cout << "Prix swap courbe : " << floatLeg - strike * fixedLeg << std::endl ;
//
//swaption
//	double mean_x_swaption	= A0 * somme_xi_swaption / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
//	double mean_x2_swaption	= A0 * A0 * somme_xi2_swaption / nbSimus; 
//	double variance_swaption = mean_x2_swaption - mean_x_swaption * mean_x_swaption ;
//
//	std::cout <<"variance : " << variance_swaption << std::endl ;
//
//	double IC_left_swaption	= mean_x_swaption - 2.57*std::sqrt(variance_swaption / nbSimus);
//	double IC_right_swaption = mean_x_swaption + 2.57*std::sqrt(variance_swaption / nbSimus);
//
//	std::cout   << "SWAPTION " << std::endl;
//	std::cout   << "MC swaption : " << mean_x_swaption << std::endl;
//	std::cout	<< "nb simulations : " << nbSimus << std::endl;
//	std::cout   << "99% confidence interval  [" << IC_left_swaption << " , " << IC_right_swaption	<< "]" << std::endl;
//
//}



//approx : show / print

//
//std::vector<double> MCforward(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer)
//{
//	double a = pApproxPricer->get_buffer_UnderlyingSwap_().get_StartDate() ;
//	double b = pApproxPricer->get_buffer_UnderlyingSwap_().get_EndDate() ;		
//param MC
//	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
//	Tenor floatTenor = pSwaption->getUnderlyingSwap().get_floatingLegTenorType() ;
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, static_cast<int>(a+b+1)) );
//	
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//
//	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
//
//	size_t fwdProbaT = static_cast<size_t>(b) ;   //end date du swap 	
//	size_t discretizationBetweenDates = 200 ;
//	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(new MC_CheyetteDD_VanillaSwaptionPricer(pModel,
//															rnGenerator, simulationStructure, fwdProbaT,
//															discretizationBetweenDates   ) ) ;
//	
//	std::vector<double> resMC = mc->price(pSwaption, nbSimus) ;
//
//	return resMC ;
//}
//
//	res[0] = prixSwaptionsPlusieursStrikesMC_QT ;
//	res[1] = IC_left ; 
//	res[2] = IC_right ; 
//	res[3] = strikes ; 
//	res[4] = moneyness ;
//std::vector<std::vector<double>> MCforward_MultipleStrikes(size_t nbSimus, 
//															CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, 
//															double	sigma_ATM)
//{
//	double a = pApproxPricer->get_buffer_UnderlyingSwap_().get_StartDate() ;
//	double b = pApproxPricer->get_buffer_UnderlyingSwap_().get_EndDate() ;		
//param MC
//	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
//	Tenor floatTenor = pSwaption->getUnderlyingSwap().get_floatingLegTenorType() ;
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, static_cast<int>(a+b+1)) );
//	
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//
//	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
//
//	size_t fwdProbaT = static_cast<size_t>(b) ;   //end date du swap 	
//	size_t discretizationBetweenDates = 200 ;
//	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(new MC_CheyetteDD_VanillaSwaptionPricer(pModel,
//															rnGenerator, simulationStructure, fwdProbaT,
//															discretizationBetweenDates   ) ) ;
//	
//pricing swaption pour plusieurs strikes 
//	double S0 = pApproxPricer->swapRate0() ; 
//	
//	std::vector<std::vector<double>> resMC = mc->priceMultipleStrikes(pSwaption, nbSimus, S0, sigma_ATM) ;
//
//	return resMC ;
//}
//
//
//void print_MCforward(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o)
//{
//
//	std::vector<double> resMC = MCforward(nbSimus, pApproxPricer) ;
//
//	double mean_x	= resMC[0] ;
//	double IC_left	= resMC[1] ;
//	double IC_right = resMC[2] ;
//
//	o   << "Monte Carlo sous proba forward Q T ; T = ;" << fwdProbaT << std::endl;
//	o   << "MC swaption : ;" << mean_x << std::endl;
//	o	<< "nb simulations : ;" << nbSimus << std::endl;
//	o   << "99% confidence interval  ;" << IC_left << " ; " << IC_right << std::endl;
//	o   << "annuity : ;" << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
//	o	<< "swap rate 0 : ;" << pApproxPricer->swapRate0() << std::endl;
//	std::cout   << "nbSimus : " << nbSimus << std::endl;
//	std::cout	<< "prix MC forward : " << mean_x << std::endl;
//	std::cout	<< "IC : [" << IC_left << " , " << IC_right << "]" << std::endl;
//	o   << nbSimus << " ; " << mean_x << " ; " << IC_left << " ; " << IC_right << std::endl;
//
//}
//
//std::vector<double> MCannuity(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer)
//{
//	
//	size_t nbPas = 200 ;
//	double a = pApproxPricer->get_buffer_T0_() ;
//	double dt = a/double(nbPas) ;
//
//	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
//	double strikeATM = pSwaption->getUnderlyingSwap().get_strike() ;
//	double S0 = pApproxPricer->get_buffer_s0_() ;
//
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	
//	std::vector<double> res(3) ;
//	double somme_xi	 = 0. ;
//	double somme_xi2 = 0. ;
//
//	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
//	{
//		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
//		std::vector<double>   gaussian_tmp(nbPas);  
//		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
//		double t = 0. ;
//		double S_t = S0 ;
//
//		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
//		{
//			double y_t, inv ;
//			if (t == 0.)
//			{
//				y_t = 0. ;
//				inv = 0. ;		
//			}
//			else 
//			{
//				y_t = pApproxPricer->get_buffer_y_bar_()(t) ; 
//				inv = pApproxPricer->inverse(t, S_t) ;
//			}
//
//			double derivee = pApproxPricer->swapRate_1stDerivative(t, inv, y_t) ;
//			double sigma = pApproxPricer->get_CheyetteDD_Model()->sigma_r(t, inv, y_t) ;
//
//			double S_t_plus_dt = S_t + derivee * sigma * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
//			S_t = S_t_plus_dt ;
//			t += dt ;
//		}
//		double res = std::max(S_t - S0, 0.) ; //  (S(T0) - K)+ 
//		
//		somme_xi += res ;
//		somme_xi2 += res * res ;
//	}
//
//	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//	double mean_x	= A0 * somme_xi / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
//	double mean_x2	= A0 * A0 * somme_xi2 / nbSimus; 
//	double variance = mean_x2 - mean_x * mean_x ;
//
//	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimus);
//	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimus);
//
//	res[0] = mean_x ;
//	res[1] = IC_left ;
//	res[2] = IC_right ;
//	 
//	return res ;
//
//}
//
//price swaptions pour différents strikes avec un seul jeu de simulations pour les differents strikes
//std::vector<std::vector<double>> MCannuityMultipleStrikes(  size_t fileNumber, size_t a, size_t b,
//															double strikeATM_Bloomberg, size_t nbSimus)
//{
//	size_t nbStrikes = 9 ;
//
//	std::vector<double> strikes(nbStrikes) ;
//
//	std::vector<double> somme_xi(nbStrikes, 0.) ;
//	std::vector<double> somme_xi2(nbStrikes, 0.) ;
//
//	CheyetteDD_Model_PTR pModel = getCalibratedModel(fileNumber, a+b) ;
//
//	std::vector<double> shiftStrike(nbStrikes) ;   //shift additif en pourcents
//	shiftStrike[0] = -2. ; shiftStrike[1] = -1. ; shiftStrike[2] = -0.5 ; shiftStrike[3] = -0.05 ;
//	shiftStrike[4] = 0. ;
//	shiftStrike[5] = 0.05 ; shiftStrike[6] = 0.5 ; shiftStrike[7] = 1. ; shiftStrike[8] = 2. ;
//
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//		strikes[i] = strike ; 
//	}
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ; 
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//	
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
//
//	size_t nbPas = 200 ;
//	double dt = a/double(nbPas) ;
//
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	
//
//	double strike = strikeATM_Bloomberg / 100. ;
//
//	VanillaSwap swap = VanillaSwap(	strike, 
//									indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
//	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer(
//					new CheyetteDD_VanillaSwaptionApproxPricer(pModel, swaption_PTR_test));
//
//	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
//	{
//		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
//		std::vector<double>   gaussian_tmp(nbPas);  
//		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
//		double t = 0. ;
//		double S_t = strikeATM_Bloomberg / 100. ;
//
//		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
//		{
//			double y_t, inv ;
//			if (t == 0.)
//			{
//				y_t = 0. ;
//				inv = 0. ;		
//			}
//			else 
//			{
//				y_t = pApproxPricer->get_buffer_y_bar_()(t) ; 
//				inv = pApproxPricer->inverse(t, S_t) ;
//			}
//
//			double derivee = pApproxPricer->swapRate_1stDerivative(t, inv, y_t) ;
//			double sigma = pApproxPricer->get_CheyetteDD_Model()->sigma_r(t, inv, y_t) ;
//
//			double S_t_plus_dt = S_t + derivee * sigma * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
//			S_t = S_t_plus_dt ;
//			t += dt ;
//		}
//
//
//		for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//		{
//			double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//
//			double res = std::max(S_t - strike, 0.) ; //  (S(T0) - K_i)+ 
//
//			somme_xi[i] += res ;
//			somme_xi2[i] += res * res ;
//		}
//
//	}
//
//	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//	double annuityA0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//	double swapRateS0 = pApproxPricer->swapRate0() ;
//
//	std::vector<double> prixSwaptionsPlusieursStrikes(nbStrikes) ;
//	std::vector<double> IC_left(nbStrikes) ;
//	std::vector<double> IC_right(nbStrikes) ;
//	std::vector<double> volBlack(nbStrikes) ;
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double mean_x = A0 * somme_xi[i] / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
//		prixSwaptionsPlusieursStrikes[i]	= mean_x ;
//		double mean_x2	= A0 * A0 * somme_xi2[i] / nbSimus; 
//		double variance = mean_x2 - mean_x * mean_x ;
//
//		IC_left[i]	= mean_x - 2.57*std::sqrt(variance / nbSimus);
//		IC_right[i] = mean_x + 2.57*std::sqrt(variance / nbSimus);
//
//
//		double volImp = NumericalMethods::Black_SwaptionImpliedVolatility(mean_x, annuityA0, swapRateS0, strike, a) ;
//		
//		volBlack[i] = volImp ;
//
//	}	
//
//	std::vector<std::vector<double>> res(5) ;
//	res[0] = strikes ;
//	res[1] = prixSwaptionsPlusieursStrikes ; 
//	res[2] = IC_left ; 
//	res[3] = IC_right ; 
//	res[4] = volBlack ; 
//
//	ecriture dans fichier
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_output_path() ;
//	fileName_s << directory << "res_smile_annuity.csv" ; 
//	std::string fileName = fileName_s.str();
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;	
//	o	<<	endl;
//
//	for (size_t i = 0 ; i < strikes.size() ; ++i)
//	{
//		o << "strike ; " << res[0][i] << "; prix " << res[1][i] << " ; " 
//															<< res[2][i] << "; " 
//															<< res[3][i] << " ; "
//															<< res[4][i] << std::endl ; 
//	}
//	
//	return res ;
//
//}
//
//
//price swaptions pour différents strikes avec un seul jeu de simulations pour les differents strikes
//std::vector<std::vector<double>> MCannuity_2_MultipleStrikes(  size_t fileNumber, size_t a, size_t b,
//															double strikeATM_Bloomberg, size_t nbSimus)
//{
//	size_t nbStrikes = 9 ;
//
//	std::vector<double> strikes(nbStrikes) ;
//
//	std::vector<double> somme_xi(nbStrikes, 0.) ;
//	std::vector<double> somme_xi2(nbStrikes, 0.) ;
//
//	CheyetteDD_Model_PTR pModel = getCalibratedModel(fileNumber, a+b) ;
//
//	std::vector<double> shiftStrike(nbStrikes) ;   //shift additif en pourcents
//	shiftStrike[0] = -2. ; shiftStrike[1] = -1. ; shiftStrike[2] = -0.5 ; shiftStrike[3] = -0.05 ;
//	shiftStrike[4] = 0. ;
//	shiftStrike[5] = 0.05 ; shiftStrike[6] = 0.5 ; shiftStrike[7] = 1. ; shiftStrike[8] = 2. ;
//
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//		strikes[i] = strike ; 
//	}
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ; 
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//	double fwdProbaT = a + b ;  
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
//
//	size_t nbPas = 200 ;
//	double dt = a/double(nbPas) ;
//
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	
//
//	double strike = strikeATM_Bloomberg / 100. ;
//
//	VanillaSwap swap = VanillaSwap(	strike, 
//									indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
//	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer(
//					new CheyetteDD_VanillaSwaptionApproxPricer(pModel, swaption_PTR_test));
//
//	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
//	{
//		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
//		std::vector<double>   gaussian_tmp(nbPas);  
//		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
//		double t = 0. ;
//		double S_t = strikeATM_Bloomberg / 100. ;
//
//		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
//		{
//			double S0 = pApproxPricer->get_buffer_s0_() ;
//			double inv = pApproxPricer->inverse(t, S0) ;
//
//			double Phi_t_s0 = pApproxPricer->Phi(t, inv) ;
//			double DPhi_t_s0 = pApproxPricer->DPhi(t, inv) ;
//
//			double approx_lineaire = Phi_t_s0 + DPhi_t_s0 * (S_t - S0) ;
//
//			double S_t_plus_dt = S_t + approx_lineaire * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
//			S_t = S_t_plus_dt ;
//			t += dt ;
//		}
//
//
//		for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//		{
//			double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//
//			double res = std::max(S_t - strike, 0.) ; //  (S(T0) - K_i)+ 
//
//			somme_xi[i] += res ;
//			somme_xi2[i] += res * res ;
//		}
//
//	}
//
//	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//	double annuityA0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//	double swapRateS0 = pApproxPricer->swapRate0() ;
//
//	std::vector<double> prixSwaptionsPlusieursStrikes(nbStrikes) ;
//	std::vector<double> IC_left(nbStrikes) ;
//	std::vector<double> IC_right(nbStrikes) ;
//	std::vector<double> volBlack(nbStrikes) ;
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double mean_x = A0 * somme_xi[i] / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
//		prixSwaptionsPlusieursStrikes[i]	= mean_x ;
//		double mean_x2	= A0 * A0 * somme_xi2[i] / nbSimus; 
//		double variance = mean_x2 - mean_x * mean_x ;
//
//		IC_left[i]	= mean_x - 2.57*std::sqrt(variance / nbSimus);
//		IC_right[i] = mean_x + 2.57*std::sqrt(variance / nbSimus);
//
//
//		double volImp = NumericalMethods::Black_SwaptionImpliedVolatility(mean_x, annuityA0, swapRateS0, strike, a) ;
//		
//		volBlack[i] = volImp ;
//
//	}	
//
//	std::vector<std::vector<double>> res(5) ;
//	res[0] = strikes ;
//	res[1] = prixSwaptionsPlusieursStrikes ; 
//	res[2] = IC_left ; 
//	res[3] = IC_right ; 
//	res[4] = volBlack ; 
//
//	ecriture dans fichier
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_output_path() ;
//	fileName_s << directory << "res_smile_annuity_2.csv" ; 
//	std::string fileName = fileName_s.str();
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;	
//	o	<<	endl;
//
//	for (size_t i = 0 ; i < strikes.size() ; ++i)
//	{
//		o << "strike ; " << res[0][i] << "; prix " << res[1][i] << " ; " 
//															<< res[2][i] << "; " 
//															<< res[3][i] << " ; "
//															<< res[4][i] << std::endl ; 
//	}
//	
//	return res ;
//
//}
//
//void print_MCannuity(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o)
//{
//
//	std::vector<double> resMCannuity = MCannuity(nbSimus, pApproxPricer) ;
//	double mean_x = resMCannuity[0] ;
//	double IC_left = resMCannuity[1] ;
//	double IC_right = resMCannuity[2] ;
//
//		o   << "Monte Carlo sous proba annuity" << std::endl;
//	o   << "MC swaption : ;" << mean_x << std::endl;
//	o	<< "nb simulations : ;" << nbSimus << std::endl;
//	o   << "99% confidence interval ; " << IC_left << " ; " << IC_right << std::endl;
//	o   << "annuity : ;" << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
//	o	<< "swap rate 0 : ;" << pApproxPricer->swapRate0() << std::endl;
//
//	std::cout   << "annuity : " << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
//	std::cout	<< "swap rate 0 : " << pApproxPricer->swapRate0() << std::endl;
//
//	std::cout   << "MC swaption : " << mean_x << std::endl;
//	std::cout	<< "nb simulations : " << nbSimus << std::endl;
//	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;
//
//	o   << nbSimus << " ; " << mean_x << " ; " << IC_left << " ; " << IC_right << std::endl;
//}
//
//
//std::vector<double> MCannuity2(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer)
//{
//	size_t nbPas = 200 ;
//	double a = pApproxPricer->get_buffer_T0_() ;
//	double dt = a/double(nbPas) ;
//
//	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
//	double strikeATM = pSwaption->getUnderlyingSwap().get_strike() ;
//
//	double S0 = pApproxPricer->swapRate0() ;		//strike et S(0)
//
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	
//	std::vector<double> res(3) ;
//	double somme_xi	 = 0. ;
//	double somme_xi2 = 0. ;
//
//	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
//	{
//		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
//		std::vector<double>   gaussian_tmp(nbPas);  
//		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
//		double t = 0. ;
//		double S_t = S0 ;
//
//		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
//		{
//			double S0 = pApproxPricer->get_buffer_s0_() ;
//			double inv = pApproxPricer->inverse(t, S0) ;
//
//			double Phi_t_s0 = pApproxPricer->Phi(t, inv) ;
//			double DPhi_t_s0 = pApproxPricer->DPhi(t, inv) ;
//
//			double approx_lineaire = Phi_t_s0 + DPhi_t_s0 * (S_t - S0) ;
//
//			double S_t_plus_dt = S_t + approx_lineaire * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
//			S_t = S_t_plus_dt ;
//			t += dt ;
//		}
//		double res = std::max(S_t - S0, 0.) ; //  (S(T0) - K)+ 
//		
//		somme_xi += res ;
//		somme_xi2 += res * res ;
//	}
//
//	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//	double mean_x	= A0 * somme_xi / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
//	double mean_x2	= A0 * A0 * somme_xi2 / nbSimus; 
//	double variance = mean_x2 - mean_x * mean_x ;
//
//	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimus);
//	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimus);
//
//	res[0] = mean_x ;
//	res[1] = IC_left ;
//	res[2] = IC_right ;
//	 
//	return res ;
//
//}
//
//void print_MCannuity2(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o)
//{
//
//	std::vector<double> resMCannuity = MCannuity2(nbSimus, pApproxPricer) ;
//	double mean_x = resMCannuity[0] ;
//	double IC_left = resMCannuity[1] ;
//	double IC_right = resMCannuity[2] ;
//
//	std::cout   << "annuity : " << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
//	std::cout	<< "swap rate 0 : " << pApproxPricer->swapRate0() << std::endl;
//
//	std::cout   << "MC swaption : " << mean_x << std::endl;
//	std::cout	<< "nb simulations : " << nbSimus << std::endl;
//	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;
//
//	o   << nbSimus << " ; " << mean_x << " ; " << IC_left << " ; " << IC_right << std::endl;
//}
//
//void MCannuity3(size_t a, size_t b, size_t nbSimus, 
//			   CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o)
//{
//	size_t nbPas = 200 ;
//	double dt = double(a)/double(nbPas) ;
//
//	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
//	double strikeATM = pSwaption->getUnderlyingSwap().get_strike() ;
//
//	double S0 = pApproxPricer->swapRate0() ;		//strike et S(0)
//	double k = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().k_ ;
//	double sigma = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().sigma_(0.5) ;  //sigma constante
//	double r0 = pApproxPricer->get_CheyetteDD_Model()->get_courbeInput_PTR()->get_f_0_t(0) ;
//
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	
//	std::vector<double> res(3) ;
//	double somme_xi	 = 0. ;
//	double somme_xi2 = 0. ;
//
//	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
//	{
//		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
//		std::vector<double>   gaussian_tmp(nbPas);  
//		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
//		double t = 0. ;
//		double S_t = S0 ;
//
//		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
//		{
//			double S0 = pApproxPricer->get_buffer_s0_() ;
//			double inv = pApproxPricer->inverse(t, S0) ;
//
//			double Phi_t_s0 = pApproxPricer->Phi(t, inv) ;
//			double DPhi_t_s0 = pApproxPricer->DPhi(t, inv) ;
//
//			double approx_lineaire = Phi_t_s0 + DPhi_t_s0 * (S_t - S0) ;
//
//			double S_t_plus_dt = S_t + approx_lineaire * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
//			S_t = S_t_plus_dt ;
//			t += dt ;
//		}
//		double res = std::max(S_t - S0, 0.) ; //  (S(T0) - K)+ 
//		
//		somme_xi += res ;
//		somme_xi2 += res * res ;
//	}
//
//	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//	double mean_x	= A0 * somme_xi / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
//	double mean_x2	= A0 * A0 * somme_xi2 / nbSimus; 
//	double variance = mean_x2 - mean_x * mean_x ;
//
//	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimus);
//	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimus);
//
//	res[0] = mean_x ;
//	res[1] = IC_left ;
//	res[2] = IC_right ;
//	 
//	o   << "Monte Carlo sous proba annuity" << std::endl;
//	o   << "MC swaption : ;" << mean_x << std::endl;
//	o	<< "nb simulations : ;" << nbSimus << std::endl;
//	o   << "99% confidence interval ; " << IC_left << " ; " << IC_right << std::endl;
//	o   << "annuity : ;" << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
//	o	<< "swap rate 0 : ;" << pApproxPricer->swapRate0() << std::endl;
//
//	std::cout   << "annuity : " << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
//	std::cout	<< "swap rate 0 : " << pApproxPricer->swapRate0() << std::endl;
//
//	std::cout   << "MC swaption : " << mean_x << std::endl;
//	std::cout	<< "nb simulations : " << nbSimus << std::endl;
//	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;
//
//	o   << nbSimus << " ; " << mean_x << " ; " << IC_left << " ; " << IC_right << std::endl;
//}
//
//
//
//
//swaption aYbY, coterminal = a + b
//calibre sur toutes les swaptions coterminales
//creation modele, et test MC vs approx
//void test_print_resultats(size_t coterminal, double k, double sigma, double m)
//{
//
////parametres
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ;  
//
//	std::vector<size_t> vectSimus(1) ; vectSimus[0] = 60000 ; //vectSimus[1] = 20000 ; vectSimus[2] = 100000 ; 
//
////ecriture dans fichier
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_output_path() ;
//	fileName_s << directory << "Comparaison_Approx_MC.csv" ; 
//	std::string fileName = fileName_s.str();
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;	
//
//	size_t fileNumber = 1 ;   //2014/01/01
//	CheyetteDD_Model_PTR pModel = getCalibratedModel(fileNumber, coterminal) ;
//
//	//CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = 
//	//						creeModeleATM(1, coterminal - 1, floatTenor, fixedTenor, k, sigma, m) ;
//	//CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
//	//pModel->print(o) ;
//
//	for (size_t a = 1 ; a < coterminal ; ++a)
//	{
//		o	<<	endl;
//		o	<<	endl;
//		o	<<	"-----------  SWAPTION " << a << "Y" << coterminal - a << "Y  -----------" << endl;
//		o	<<	endl;
//
//
//		Tenor floatTenor = Tenor::_6M ;
//		Tenor fixedTenor = Tenor::_1YR ; 
//		double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//		
//		size_t indexStart = size_t(a / tenor) ;
//		size_t indexEnd = size_t(coterminal / tenor) ;
//		double fwdProbaT = coterminal ;  
//		LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, coterminal + 1) );
//
//		double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//
//		VanillaSwap swap = VanillaSwap(	strike, 
//										indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
//		VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//
//		CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer(
//						new CheyetteDD_VanillaSwaptionApproxPricer(pModel, swaption_PTR_test));
//
//		CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = 
//								creeModeleATM(a, coterminal - a, floatTenor, fixedTenor, k, sigma, m) ;
//
//		//approx
//		o   << "prix par approximation ;" << pApproxPricer->prixSwaptionApproxPiterbarg() << std::endl;
//
//		////MC forward QT
//		//o	<<	endl;
//		//o   << "Monte Carlo sous proba forward Q T ; T = ;" << coterminal << std::endl;
//		//o   << "nbSimus ; MC swaption ; IC_left ; IC_right" << std::endl;
//		//for (size_t indexSimu = 0 ; indexSimu < vectSimus.size() ; ++indexSimu)
//		//{
//		//	MCforward(vectSimus[indexSimu], pApproxPricer) ;	
//		//}
//
//		////MC annuity
//		//o	<<	endl;
//		//o   << "Monte Carlo sous proba annuity" << std::endl;
//		//o   << "nbSimus ; MC swaption ; IC_left ; IC_right" << std::endl;
//		//for (size_t indexSimu = 0 ; indexSimu < vectSimus.size() ; ++indexSimu)
//		//{
//		//	print_MCannuity(vectSimus[indexSimu], pApproxPricer, o) ; 
//		//}
//		o	<<	endl;
//
//	}
//
//	o.close() ;
//}
//
//
//evaluation de Phi(t, S_t) aux points t = 0...T, S_t = 0, ...
//void evaluatePhi(double t_min, double t_max, double S_min, double S_max)
//{
//ecriture dans fichier
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_output_path() ;
//	fileName_s << directory << "Phi.csv" ; 
//	std::string fileName = fileName_s.str();
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;	
//	o	<<	endl;
//
//	size_t fileNumber = 24 ;   //2014/01/01
//	size_t coterminal = 10 ;
//	CheyetteDD_Model_PTR pModel = getCalibratedModel(fileNumber, coterminal) ;
//
//	double dt = 0.01 ;
//	double dS = 0.001 ;
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ;  
//	double k = 0.02 ;
//	double sigma = 0.8 ;
//	double m = 0.15 ;
//	size_t a = 1 ;
//	size_t b = 11 ;
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = 
//							creeModeleATM(a, b, floatTenor, fixedTenor, k, sigma, m) ;
//
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = creeModeleATM(a, b, 
//														 pModel,
//														 floatTenor, fixedTenor) ;
//
//	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
//	pModel->print(o) ;
//
//	size_t nb_dt = static_cast<size_t>((t_max - t_min) / dt) ;
//	size_t nb_dS = static_cast<size_t>((S_max - S_min) / dS) ;
//	o	<< " ; ; t : " ;
//	for (size_t i = 0 ; i <= nb_dt ; ++i)
//	{
//		o << t_min + i * dt << " ; "  ;		
//	}
//	o << endl;
//	for (size_t indexS = 0 ; indexS <= nb_dt ; ++indexS)
//	{
//		o	<< "S :  ; " << S_min + indexS * dS << " ; " ;
//		for (size_t indexT = 0 ; indexT <= nb_dt ; ++indexT)
//		{
//			double t = t_min + indexT * dt ;
//			double S_t = S_min + indexS * dS ;
//
//			double	y_t = pApproxPricer->get_buffer_y_bar_()(t) ;
//			double	inv = pApproxPricer->inverse(t, S_t) ;
//			
//			double derivee = pApproxPricer->swapRate_1stDerivative(t, inv, y_t) ;
//			double sigma = pApproxPricer->get_CheyetteDD_Model()->sigma_r(t, inv, y_t) ;
//
//			std::cout << "t : " << t << ", S_t : " << S_t << std::endl ;
//			o	<< derivee * sigma << " ; " ;		
//		}
//		o << endl;
//	}
//	std::cout << "OK" << std::endl ;
//	o.close() ;
//}
//
//
//void evaluatePhi_lineaire(double t_min, double t_max, double S_min, double S_max)
//{
//ecriture dans fichier
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_output_path() ;
//	fileName_s << directory << "Phi_lineaire.csv" ; 
//	std::string fileName = fileName_s.str();
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;	
//	o	<<	endl;
//
//	size_t fileNumber = 24 ;   //2014/01/01
//	size_t coterminal = 10 ;
//	CheyetteDD_Model_PTR pModel = getCalibratedModel(fileNumber, coterminal) ;
//
//	double dt = 0.01 ;
//	double dS = 0.001 ;
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ;  
//	double k = 0.02 ;
//	double sigma = 0.8 ;
//	double m = 0.15 ;
//	size_t a = 1 ;
//	size_t b = 11 ;
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = 
//							creeModeleATM(a, b, floatTenor, fixedTenor, k, sigma, m) ;
//
//	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = creeModeleATM(a, b, 
//														 pModel,
//														 floatTenor, fixedTenor) ;
//
//
//
//
//	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
//	pModel->print(o) ;
//
//	size_t nb_dt = static_cast<size_t>((t_max - t_min) / dt) ;
//	size_t nb_dS = static_cast<size_t>((S_max - S_min) / dS) ;
//	o	<< " ; ; t : " ;
//	for (size_t i = 0 ; i <= nb_dt ; ++i)
//	{
//		o << t_min + i * dt << " ; "  ;		
//	}
//	o << endl;
//	for (size_t indexS = 0 ; indexS <= nb_dt ; ++indexS)
//	{
//		o	<< "S :  ; " << S_min + indexS * dS << " ; " ;
//		for (size_t indexT = 0 ; indexT <= nb_dt ; ++indexT)
//		{
//			double t = t_min + indexT * dt ;
//			double S_t = S_min + indexS * dS ;
//
//			double s_bar = pApproxPricer->get_buffer_s0_() ;
//
//			double	inv = pApproxPricer->inverse(t, s_bar) ;
//
//			double phi_t_s_bar = pApproxPricer->Phi(t, inv) ;
//			double Dphi_t_s_bar = pApproxPricer->DPhi(t, inv) ;
//
//			std::cout << "t : " << t << ", S_t : " << S_t << std::endl ;
//			o	<< phi_t_s_bar + Dphi_t_s_bar * (S_t - s_bar) << " ; " ;		
//		}
//		o << endl;
//	}
//	std::cout << "OK" << std::endl ;
//	o.close() ;
//}
//
//swaption aY bY... 
//void smile(size_t a, size_t b, double strikeATM_Bloomberg)
//{
//calibre sur le modèle désiré et retourne le modele calibre
//
//	size_t fileNumber = 1 ;   //2008/04/02
//	size_t coterminal = a + b ; 
//	CheyetteDD_Model_PTR pModel = getCalibratedModel(fileNumber, coterminal) ;
//
//	std::vector<double> shiftStrike(9) ;   //shift additif en pourcents
//	shiftStrike[0] = -2. ; shiftStrike[1] = -1. ; shiftStrike[2] = -0.5 ; shiftStrike[3] = -0.05 ;
//	shiftStrike[4] = 0. ;
//	shiftStrike[5] = 0.05 ; shiftStrike[6] = 0.5 ; shiftStrike[7] = 1. ; shiftStrike[8] = 2. ;
//
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_output_path() ;
//	fileName_s << directory << "smile_new.csv" ; 
//	std::string fileName = fileName_s.str();
//
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;
//
//	pModel->print(o) ;
//	o << endl ;
//	pModel->get_courbeInput_PTR()->printHorizontal(o) ;
//	o << endl ;
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ; 
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//	double fwdProbaT = b ;  //size_t vers double ok
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
//	
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		o << strikeATM_Bloomberg + shiftStrike[i] << " ; " ; 
//	}
//	o << endl ;
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//
//		VanillaSwap swap = VanillaSwap(	strike, 
//										indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
//		VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//
//		CheyetteDD_VanillaSwaptionApproxPricer approx(pModel, swaption_PTR_test);
//
//		double approxPrice = approx.prixSwaptionApproxPiterbarg() ;
//
//		double annuityA0 = approx.swapRateDenominator(0., 0., 0.) ;
//		double swapRateS0 = approx.swapRate0() ;
//		double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 															
//							swapRateS0, strike, a) ;
//				
//		o << volBlack << ";" ;
//	}
//	o << endl ;
//	o.close();
//						 
//}
//
//
//void helpPrinter(const std::string& descriptionData, const std::vector<double>& data, std::ofstream& o)
//{
//	o << descriptionData << " ; " ;
//	for (size_t i = 0 ; i < data.size() ; ++i)
//	{
//		o << data[i] << " ; " ; 
//	}
//	o << endl ;
//}
//
//
//void print_smile(size_t a, size_t b, double strikeATM_Bloomberg, 
//				 double price, double annuity0, double swapRate0, std::ofstream& o)
//{
//	std::vector<double> shiftStrike(9) ;   //shift additif en pourcents
//	shiftStrike[0] = -2./100. ; shiftStrike[1] = -1./100. ; shiftStrike[2] = -0.5/100. ; shiftStrike[3] = -0.05/100. ;
//	shiftStrike[4] = 0. ;
//	shiftStrike[5] = 0.05/100. ; shiftStrike[6] = 0.5/100. ; shiftStrike[7] = 1./100. ; shiftStrike[8] = 2./100. ;
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ; 
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//	double fwdProbaT = b ;  //size_t vers double ok
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
//	
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		o << strikeATM_Bloomberg + shiftStrike[i] << " ; " ; 
//	}
//	o << endl ;
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//		double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(price, annuity0, 															
//							swapRate0, strike, a) ;
//				
//		o << volBlack << ";" ;
//	}
//	o << endl ;
//	o.close();			 //flux appele dans les fichiers suivants
//}
//
//with standardized moneyness
//void print_smile_standard_moneyness(size_t a, size_t b, double strikeATM_Bloomberg, 
//									double price, double annuity0, double swapRate0, std::ofstream& o)
//{
//	std::vector<double> shiftStrike(9) ;   //shift additif en pourcents
//	shiftStrike[0] = -2./100. ; shiftStrike[1] = -1./100. ; shiftStrike[2] = -0.5/100. ; shiftStrike[3] = -0.05/100. ;
//	shiftStrike[4] = 0. ;
//	shiftStrike[5] = 0.05/100. ; shiftStrike[6] = 0.5/100. ; shiftStrike[7] = 1./100. ; shiftStrike[8] = 2./100. ;
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ; 
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//	double fwdProbaT = b ;  //size_t vers double ok
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
//	
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		o << strikeATM_Bloomberg + shiftStrike[i] << " ; " ; 
//	}
//	o << endl ;
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//		double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(price, annuity0, 															
//							swapRate0, strike, a) ;
//				
//		o << volBlack << ";" ;
//	}
//	o << endl ;
//	o.close();			 //flux appele dans les fichiers suivants
//}
//1er element du vecteur : strikes
//2eme element du vecteur : prix
//3eme element du vecteur : vol Black
//pour un strike de 2%, rentrer 2
//std::vector<std::vector<double>> smileVect(size_t fileNumber, size_t a, size_t b,
//												  double strikeATM_Bloomberg, size_t nbSimus)
//{
//
//	size_t nbStrikes = 9 ;
//
//	std::vector<double> strikes(nbStrikes) ;
//	std::vector<double> prixSwaptionsPlusieursStrikes(nbStrikes) ;
//	std::vector<double> volBlack(nbStrikes) ;
//
//	CheyetteDD_Model_PTR pModel = getCalibratedModel(fileNumber, a+b) ;
//
//	std::vector<double> shiftStrike(nbStrikes) ;   //shift additif en pourcents
//	shiftStrike[0] = -2. ; shiftStrike[1] = -1. ; shiftStrike[2] = -0.5 ; shiftStrike[3] = -0.05 ;
//	shiftStrike[4] = 0. ;
//	shiftStrike[5] = 0.05 ; shiftStrike[6] = 0.5 ; shiftStrike[7] = 1. ; shiftStrike[8] = 2. ;
//
//	Tenor floatTenor = Tenor::_6M ;
//	Tenor fixedTenor = Tenor::_1YR ; 
//	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
//	
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//	double fwdProbaT = b ;  
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
//	
//	for (size_t i = 0 ; i < shiftStrike.size() ; ++i)
//	{
//		double strike = (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//		strikes[i] = strike ; 
//
//		VanillaSwap swap = VanillaSwap(	strike, indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
//		VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//
//		CheyetteDD_VanillaSwaptionApproxPricer_PTR approx(
//			new CheyetteDD_VanillaSwaptionApproxPricer( pModel, swaption_PTR_test)) ;
//
//		double MCprice = MCannuity(nbSimus, approx)[0] ;
//		
//		prixSwaptionsPlusieursStrikes[i] = MCprice ;
//
//		double annuityA0 = approx->swapRateDenominator(0., 0., 0.) ;
//		double swapRateS0 = approx->swapRate0() ;
//		double volImp = NumericalMethods::Black_SwaptionImpliedVolatility(MCprice, annuityA0, swapRateS0, strike, a) ;
//		
//		volBlack[i] = volImp ;
//	}
//
//	std::vector<std::vector<double>> res(3) ;
//	res[0] = strikes ;
//	res[1] = prixSwaptionsPlusieursStrikes ; 
//	res[2] = volBlack ; 
//
//	for (size_t i = 0 ; i < strikes.size() ; ++i)
//	{
//		std::cout << "strike : " << res[0][i] << " , prix " << res[1][i] << " , " << res[2][i] << std::endl ; 
//	}
//	
//	return res ;
//}
//
//print matrix 
//
//void test_printMatrix(std::vector<size_t> expiry_maturity, double strike, Tenor floatingLegTenor, Tenor fixedLegTenor, 
//					 int curveChoice, int shiftChoice, 
//					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
//					 std::vector<size_t> nbSimus)
//{	
////courbe spot
//	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
//	
//	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
//	Piecewiseconst_RR_Function m(x, m_y) ;
//	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
//
//	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
//	modele_test_PTR->show() ;
//
////swaption grid
//	size_t dimMatrix = expiry_maturity.size() ;
//	Matrice matriceApprox (dimMatrix, dimMatrix) ;
//	Matrice matriceS0 (dimMatrix, dimMatrix) ;
//	Matrice matriceBbarre (dimMatrix, dimMatrix) ;
//	Matrice matriceVolBlack (dimMatrix, dimMatrix) ;
//	Matrice matriceMCprice (dimMatrix, dimMatrix) ;
//	Matrice matriceMC_ICinf (dimMatrix, dimMatrix) ;
//	Matrice matriceMC_ICsup (dimMatrix, dimMatrix) ;
//	Matrice matriceRelativeError (dimMatrix, dimMatrix) ;
//	Matrice matriceAbsoluteError (dimMatrix, dimMatrix) ;
//	Matrice matriceKtilde (dimMatrix, dimMatrix) ;
//	for (size_t index1 = 0 ; index1 < expiry_maturity.size() ; ++index1)
//	{
//		for (size_t index2 = 0 ; index2 < expiry_maturity.size() ; ++index2)
//		{
//			test_printElement(strike, index1, index2, expiry_maturity, floatingLegTenor, fixedLegTenor, curveChoice, shiftChoice, 
//									modele_test_PTR, nbSimus,
//									matriceApprox, matriceS0, matriceBbarre, matriceVolBlack, 
//									matriceMCprice, matriceMC_ICinf, matriceMC_ICsup, 
//									matriceRelativeError, matriceAbsoluteError, matriceKtilde) ;	
//		}	
//		std::cout << "index 1 : " << index1 << std::endl ;
//	}
//
//
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_runtime_datapath() ;
//	fileName_s << directory << "test.csv" ; 
//	std::string fileName = fileName_s.str();
//
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;
//	modele_test_PTR->print(o) ;
//	o	<<	endl;
//	o	<< "strike choisi : " << strike << endl ;
//	switch (shiftChoice)
//	{
//		case 1 :
//			o << "shift 1 : r(t) / r(0)" << endl ;
//			break ;
//		case 2 :
//			o << "shift 2 : r(t) / f(0, t)" << endl ;
//			break ;
//		case 3 :
//			o << "shift 3 : S(t) / S(0)" << endl ;
//			break ;
//		default:
//			std::cout << "invalide shiftChoice" << std::endl ;
//			throw "exception" ;
//			break;
//	}
//	
//	courbe_PTR_test->print(o) ;
//	o << endl ;
//
//	matrixPrint(o, expiry_maturity, "approximation price", matriceApprox) ;
//	matrixPrint(o, expiry_maturity, "MC price - nbSimu ", matriceMCprice) ;
//	matrixPrint(o, expiry_maturity, "Erreur relative", matriceRelativeError) ;
//	matrixPrint(o, expiry_maturity, "Erreur absolue : MCprice - approx", matriceAbsoluteError) ;
//	matrixPrint(o, expiry_maturity, "S0 (strike atm)", matriceS0) ;
//	matrixPrint(o, expiry_maturity, "b barre", matriceBbarre) ;
//	matrixPrint(o, expiry_maturity, "vol Black", matriceVolBlack) ;
//	matrixPrint(o, expiry_maturity, "MC IC inf", matriceMC_ICinf) ;
//	matrixPrint(o, expiry_maturity, "MC IC sup", matriceMC_ICsup) ;
//	matrixPrint(o, expiry_maturity, "K tilde", matriceKtilde) ;
//
//	o.close();
//}
//
//void test_printElement(double strike, size_t index1, size_t index2, std::vector<size_t> expiry_maturity, 
//						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
//						 int curveChoice, int shiftChoice, 
//						 CheyetteDD_Model_PTR modele_test_PTR, 
//						 std::vector<size_t> nbSimus, 
//						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
//						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
//						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde)
//{	
//	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
//	
//	size_t a = expiry_maturity[index1] ;
//	size_t b = expiry_maturity[index2] ;
//
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//	size_t fwdProbaT = a + b ; 
//
////param MC
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
//	size_t	discretizationBetweenDates = 200 ;
//
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(
//					new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR, rnGenerator,
//															simulationStructure, fwdProbaT,
//															discretizationBetweenDates   ) ) ;
//
//	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
//	
//	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);
//
//	std::cout << "prixSwaption" << std::endl ;
//	double approxPrice = approx.prixSwaptionApproxPiterbarg() ;
//	double b_barre = approx.get_buffer_b_barre_() ;
//	std::cout << "approxPrice : " << approxPrice << std::endl ;
//	std::cout << "b barre :     " << b_barre << std::endl ;
//	std::cout << "  " << std::endl ;
//	matriceApprox (index1, index2)  = approxPrice ; 
//	matriceBbarre (index1, index2)  = b_barre ; 
//
////MC
//	//double nominal = 1 ;
//	//LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, a+b+1)) ;
//	//GeneticSwap_CONSTPTR genericSwapTest = InstrumentFactory::createVanillaSwap(strikeATM, indexStart, 
//	//												indexEnd, floatingLegTenor, fixedLegTenor,
//	//												swapStructure, nominal);
//
//	//size_t maturity = indexStart ; //index
//	//GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;
//
//	size_t nbMC = nbSimus.size() ;
//	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
//	for (size_t i = 0 ; i < nbMC ; ++i)
//	{
//		//std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;
//		std::vector<double> resMC = mc->price(swaption_PTR_test, nbSimus[i]) ;
//		matriceMCprice (index1, index2)  = resMC[0] ;
//		matriceMC_ICinf (index1, index2) = resMC[1] ;
//		matriceMC_ICsup (index1, index2) = resMC[2] ;
//	}
//
//	double annuityA0 = approx.swapRateDenominator(0., 0., 0.) ;
//	double swapRateS0 = approx.swapRate0() ;
//	matriceS0(index1, index2) = swapRateS0 ;
//	matriceKtilde(index1, index2) = b_barre * strike + (1 - b_barre) * swapRateS0 ;
//	double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
//																	swapRateS0, strike, a) ;
//	matriceVolBlack (index1, index2) = volBlack ;
//}
//
//
////param 
////floating leg tenor vs fixed leg tenor (pas vraiment un tenor mais frequence des flux)
////size(x) = size(m_y) + 1 = size(sigma_y) + 1
//void test_printMatrix_ATM(std::vector<size_t> expiry_maturity, Tenor floatingLegTenor, Tenor fixedLegTenor, 
//					 int curveChoice, int shiftChoice, 
//					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
//					 std::vector<size_t> nbSimus)
//{	
////courbe spot
//	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
//	
//	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
//	Piecewiseconst_RR_Function m(x, m_y) ;
//	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
//
//	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
//	modele_test_PTR->show() ;
//
////swaption grid
//	size_t dimMatrix = expiry_maturity.size() ;
//	Matrice matriceApprox (dimMatrix, dimMatrix) ;
//	Matrice matriceS0 (dimMatrix, dimMatrix) ;
//	Matrice matriceBbarre (dimMatrix, dimMatrix) ;
//	Matrice matriceVolBlack (dimMatrix, dimMatrix) ;
//	Matrice matriceMCprice (dimMatrix, dimMatrix) ;
//	Matrice matriceMC_ICinf (dimMatrix, dimMatrix) ;
//	Matrice matriceMC_ICsup (dimMatrix, dimMatrix) ;
//	Matrice matriceRelativeError (dimMatrix, dimMatrix) ;
//	Matrice matriceAbsoluteError (dimMatrix, dimMatrix) ;
//	Matrice matriceKtilde (dimMatrix, dimMatrix) ;
//	for (size_t index1 = 0 ; index1 < expiry_maturity.size() ; ++index1)
//	{
//		for (size_t index2 = 0 ; index2 < expiry_maturity.size() ; ++index2)
//		{
//			test_printElement_ATM(index1, index2, expiry_maturity, floatingLegTenor, fixedLegTenor, curveChoice, shiftChoice, 
//									modele_test_PTR, nbSimus,
//									matriceApprox, matriceS0, matriceBbarre, matriceVolBlack, 
//									matriceMCprice, matriceMC_ICinf, matriceMC_ICsup, 
//									matriceRelativeError, matriceAbsoluteError, matriceKtilde) ;	
//		}	
//		std::cout << "index 1 : " << index1 << std::endl ;
//	}
//
//
//	std::stringstream fileName_s ;
//	std::string directory = LMMPATH::get_runtime_datapath() ;
//	fileName_s << directory << "test_shift.csv" ; 
//	std::string fileName = fileName_s.str();
//
//	ofstream o;
//	o.open(fileName,  ios::out | ios::app );
//	o	<<	endl;
//	modele_test_PTR->print(o) ;
//	o	<<	endl;
//	switch (shiftChoice)
//	{
//		case 1 :
//			o << "shift 1 : r(t) / r(0)" << endl ;
//			break ;
//		case 2 :
//			o << "shift 2 : r(t) / f(0, t)" << endl ;
//			break ;
//		case 3 :
//			o << "shift 3 : S(t) / S(0)" << endl ;
//			break ;
//		default:
//			std::cout << "invalide shiftChoice" << std::endl ;
//			throw "exception" ;
//			break;
//	}
//	
//	courbe_PTR_test->print(o) ;
//	o << endl ;
//	
//
//	matrixPrint(o, expiry_maturity, "approximation price ATM", matriceApprox) ;
//	matrixPrint(o, expiry_maturity, "MC price - nbSimu ", matriceMCprice) ;
//	matrixPrint(o, expiry_maturity, "Erreur relative", matriceRelativeError) ;
//	matrixPrint(o, expiry_maturity, "Erreur absolue : MCprice - approx", matriceAbsoluteError) ;
//	matrixPrint(o, expiry_maturity, "S0 (strike atm)", matriceS0) ;
//	matrixPrint(o, expiry_maturity, "b barre", matriceBbarre) ;
//	matrixPrint(o, expiry_maturity, "vol Black", matriceVolBlack) ;
//	matrixPrint(o, expiry_maturity, "MC IC inf", matriceMC_ICinf) ;
//	matrixPrint(o, expiry_maturity, "MC IC sup", matriceMC_ICsup) ;
//	matrixPrint(o, expiry_maturity, "K tilde", matriceKtilde) ;
//
//	o.close();
//}
//
//void matrixPrint(std::ofstream& o, std::vector<size_t> expiry_maturity, std::string chaine, const Matrice& M)
//{
//	o << endl ;
//	o << endl ;
//	o << chaine << endl ;
//	o << "expiry (a) \\ tenor swap (b) ; " ;
//		for (size_t i = 0 ; i < expiry_maturity.size() ; i++)
//		{
//			o << expiry_maturity[i] << "Y ;" ;
//		}
//		o	<<	endl;
//		MatrixPrintElement<Matrice> matrixPrinter("azzrehrj", M) ;
//		for (size_t rowIndex = 0 ; rowIndex < expiry_maturity.size() ; rowIndex++)
//		{
//			o << expiry_maturity[rowIndex] << "Y ;" ;
//			matrixPrinter.print_element2(o, rowIndex) ;
//			o << endl ;
//		}
//}
//
//void test_printElement_ATM(size_t index1, size_t index2, std::vector<size_t>& expiry_maturity, 
//						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
//						 int curveChoice, int shiftChoice, 
//						 CheyetteDD_Model_PTR modele_test_PTR, 
//						 std::vector<size_t>& nbSimus, 
//						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
//						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
//						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde)
//{	
//	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
//	
//	size_t a = expiry_maturity[index1] ;
//	size_t b = expiry_maturity[index2] ;
//
//	size_t indexStart = size_t(a / tenor) ;
//	size_t indexEnd = size_t((a+b) / tenor) ;
//
////approx
//	double strike = 1000. ; //temporaire avant de mettre strike ATM
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
//	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
//	
//	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);
//
////calcul du strike ATM pour le swap
//	double strikeATM = approx.swapRate0() ;
//	swap.set_strike(strikeATM) ;
//	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
//	CheyetteDD_VanillaSwaptionApproxPricer approxATM(modele_test_PTR, swaption_PTR_test_ATM);
//
//
////param MC
//	size_t fwdProbaT = a + b ;
//	size_t	discretizationBetweenDates = 200 ;
//
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(
//					new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR, rnGenerator,
//															simulationStructure, fwdProbaT,
//															discretizationBetweenDates   ) ) ;
//
//
//	std::cout << "prixSwaption" << std::endl ;
//	double approxPrice = approxATM.prixSwaptionApproxPiterbarg() ;
//	double b_barre = approxATM.get_buffer_b_barre_() ;
//	std::cout << "approxPrice : " << approxPrice << std::endl ;
//	std::cout << "b barre :     " << b_barre << std::endl ;
//	std::cout << "  " << std::endl ;
//	matriceApprox (index1, index2)  = approxPrice ; 
//	matriceBbarre (index1, index2)  = b_barre ; 
//
////MC
//	//double nominal = 1 ;
//	//LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, a+b+1)) ;
//	//GeneticSwap_CONSTPTR genericSwapTest = InstrumentFactory::createVanillaSwap(strikeATM, indexStart, 
//	//												indexEnd, floatingLegTenor, fixedLegTenor,
//	//												swapStructure, nominal);
//
//	//size_t maturity = indexStart ; //index
//	//GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;
//
//	size_t nbMC = nbSimus.size() ;
//	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
//	for (size_t i = 0 ; i < nbMC ; ++i)
//	{
//		//std::vector<double> resMC = mc->price(genericSwaptionTest_ATM, nbSimus[i]) ;
//		std::vector<double> resMC = mc->price(swaption_PTR_test_ATM, nbSimus[i]) ;
//		vectPrixMC[i] = resMC[0] ;
//		vectICinf[i]  = resMC[1] ;
//		vectICsup[i]  = resMC[2] ;
//		matriceMCprice (index1, index2)  = resMC[0] ;
//		matriceMC_ICinf (index1, index2) = resMC[1] ;
//		matriceMC_ICsup (index1, index2) = resMC[2] ;
//	}
//
//	matriceAbsoluteError(index1, index2) = matriceMCprice(index1, index2) - matriceApprox(index1, index2) ;
//	matriceRelativeError(index1, index2) = matriceAbsoluteError(index1, index2) / matriceMCprice(index1, index2) ;
//
//	double annuityA0 = approxATM.swapRateDenominator(0., 0., 0.) ;
//	double swapRateS0 = approxATM.swapRate0() ;
//	matriceS0(index1, index2) = swapRateS0 ;
//
//	double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
//																	swapRateS0, strikeATM, a) ;
//	matriceVolBlack (index1, index2) = volBlack ;
//	matriceKtilde (index1, index2) = b_barre * strikeATM + (1 - b_barre) * swapRateS0 ; 
//	//mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
//	//						a, b, genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
//
//}
//
//
