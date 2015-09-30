#include <Cheyette\unit_test\TestApproxCall.h>



Cheyette_SwaptionPricer_QuadApprox_PTR createQuadApproxPricer_PTR(size_t xmax,  int numCourbe, 
																double k, double aValue, double bValue, double cValue)
{
	//a(t) + b(t) * x(t) + c(t) * x(t)^2

	CheyetteQuad_Model_PTR	cheyetteQuad_Model_Test_PTR(creeCheyetteQuad_Modele_PTR(xmax, numCourbe, 
																					k, aValue, bValue, cValue) ) ;

	Cheyette_SwaptionPricer_QuadApprox_PTR approxPricerTest_PTR(
							new Cheyette_SwaptionPricer_QuadApprox(	cheyetteQuad_Model_Test_PTR)); 
	return approxPricerTest_PTR ;
}

Cheyette_SwaptionPricer_QuadApprox_PTR createQuadApproxPricer_PTR(size_t xmax,  CourbeInput_PTR pCourbeInput, 
																double k, double aValue, double bValue, double cValue)
{
	//a(t) + b(t) * x(t) + c(t) * x(t)^2

	CheyetteQuad_Model_PTR	cheyetteQuad_Model_Test_PTR(creeCheyetteQuad_Modele_PTR(xmax, pCourbeInput, 
																					k, aValue, bValue, cValue) ) ;

	Cheyette_SwaptionPricer_QuadApprox_PTR approxPricerTest_PTR(
							new Cheyette_SwaptionPricer_QuadApprox(	cheyetteQuad_Model_Test_PTR)); 
	return approxPricerTest_PTR ;
}

//double phi(double St){	return St ;}
//double phiP(double St){	return 1. ;}   //phi prime
//double phiS(double St){	return 0. ;}	 //phi seconde
//double phiMoinsUn(double St){return 1. / St ;}

double phi(double St){	return sqrt(St) ;}
double phiP(double St){	return 1. / (2. * sqrt(St)) ;}   //phi prime
double phiS(double St){	return - 1./ 4. * pow(St, - 3./2.) ;}	 //phi seconde
double phiMoinsUn(double St){return 1. / sqrt(St) ;}

//approximation call


void testCallMC()
{
	size_t nbSimus = 1000000 ;		//vol de 100%, mettre beaucoup de simus >= 500 000 
	double strike = 100. ;			//prendre ATM !!
	double S0 = 100. ;
	double T = 1. ;
	std::vector<double> resMC = approxMC_call(nbSimus, strike, S0, T) ;

	//test : ... retour une variable locale
	std::cout << resMC[0] << ", " << resMC[1]  << std::endl ;
}

void testCallApprox()
{
	int numCourbe = 2 ;
	size_t xmax = 2 ;
	double aValue = 0.2 ;
	double bValue = 0.05 ;
	double cValue = 0.01 ;
	double k = 0.02 ;
	Cheyette_SwaptionPricer_QuadApprox_PTR approxPricerPTR = createQuadApproxPricer_PTR(xmax, numCourbe, 
																						k, aValue, bValue, cValue) ;

	//test des fonctions ATM (regle de l'Hopital)
	std::cout << "omega 0 ATM :  " << approxPricerPTR->omega0_ATM() << std::endl ;
	double epsilon = 0.1 ; //1e-6 ;
	std::cout << "omega 0 100 - 2e :  " << approxPricerPTR->omega0(100 - 2 * epsilon) << std::endl ;
	std::cout << "omega 0 100 - e :  " << approxPricerPTR->omega0(100 - epsilon) << std::endl ;
	std::cout << "omega 0 100 + e :  " << approxPricerPTR->omega0(100 + epsilon) << std::endl ;
	std::cout << "omega 0 100 + 2e :  " << approxPricerPTR->omega0(100 + 2 * epsilon) << std::endl ;

	std::cout << "omega 1 ATM :  " << approxPricerPTR->omega1_ATM() << std::endl ;
	std::cout << "omega 1 100 - 3e :  " << approxPricerPTR->omega1(100 - 3 * epsilon) << std::endl ;
	std::cout << "omega 1 100 - 2e :  " << approxPricerPTR->omega1(100 - 2 * epsilon) << std::endl ;
	std::cout << "omega 1 100 - e :  " << approxPricerPTR->omega1(100 - epsilon) << std::endl ;
	std::cout << "omega 1 100 + e :  " << approxPricerPTR->omega1(100 + epsilon) << std::endl ;
	std::cout << "omega 1 100 + 10e :  " << approxPricerPTR->omega1(100 + 10 * epsilon) << std::endl ;

	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;

	VanillaSwaption_PTR vanillaSwaption_Test_PTR = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;
	std::cout << "prix approx quadratique : " << approxPricerPTR->price(vanillaSwaption_Test_PTR) << std::endl ;	
}


//dS(t) = Phi(S_t) dW_t
//S(0) = S0
std::vector<double> approxMC_call(const size_t nbSimus, const double strike, const double S0, const double T) 
{
	boost::function<double(double)> f1 = boost::bind(phi, _1);
	Boost_RR_Function_PTR funcPhi (new Boost_RR_Function(f1)) ;

	std::vector<double> res(3) ;
	double somme_xi   = 0.0;
	double somme_xi_2 = 0.0;

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	for(size_t itrSimulation=0; itrSimulation < nbSimus ; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimus == 0){std::cout << double(itrSimulation)/double(nbSimus)*100 << "%" << std::endl ;}	

		//simulation de S_0 jusqu'à T
		size_t nbStepsPerYear = 250 ; 
		double dt = 1. / static_cast<double>(nbStepsPerYear) ;
		std::vector<double>   gaussian_tmp(static_cast<size_t>(T * nbStepsPerYear));  
		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.

		double t = 0. ;
		double S_t = S0 ;
		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= T * nbStepsPerYear ; ++pasDiscretisation)    
		{
			double t_plus_dt = t + dt ;
			double S_t_plus_dt =  S_t + funcPhi->operator()(S_t) * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
						
			S_t = S_t_plus_dt ;
			t = t_plus_dt ;
		}
			
		double payoff = std::max(S_t - strike, 0.) ;  
		
		somme_xi	+= payoff ;	
		somme_xi_2	+= payoff * payoff ;
	}

	double mean_x	= somme_xi / nbSimus ; 
	double mean_x2	= somme_xi_2 / nbSimus ; 
	double variance = mean_x2 - mean_x * mean_x ;

	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimus);
	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimus);

	res[0] = mean_x ;
	res[1] = IC_left ;
	res[2] = IC_right ;

	std::cout   << "prix MC swaption : " << mean_x << std::endl;
	std::cout	<< "nbSimulations    : " << nbSimus << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;

	return res ;

}

//MC contre approx quadratique
void test_MC_approx()
{
	int numCourbe = 2 ;
	size_t xmax = 3 ;
	double aValue = 0.1 ;
	double bValue = 0.05 ;
	double cValue = 0.01 ;
	double k = 0.02 ;

	CheyetteQuad_Model_PTR cheyetteQuad_Model_PTR = creeCheyetteQuad_Modele_PTR(xmax, numCourbe, 
																				k, aValue, bValue, cValue) ;
	Cheyette_SwaptionPricer_QuadApprox_PTR pQuadApproxPricer = createQuadApproxPricer_PTR(xmax, numCourbe, 
																				k, aValue, bValue, cValue) ;

	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 6; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;

	VanillaSwaption_PTR vanillaSwaption_Test_PTR = createSwaption(strike, indexStart, indexEnd, 
																  floatingLegTenorType, fixedLegTenorType) ;
	std::cout << "prix approx quadratique : " << pQuadApproxPricer->price(vanillaSwaption_Test_PTR) << std::endl ;	

	size_t fwdProbaT = 4;
	size_t discretizationBetweenDates = 200 ;
	LMMTenorStructure_PTR	pTenorStructure(new LMMTenorStructure(Tenor::_6M, 4)) ;
	MC_Cheyette_VanillaSwaptionPricer_PTR mcQuadPricer = creeMC_SwaptionPricer_PTR(cheyetteQuad_Model_PTR, 
																	pTenorStructure,
																	fwdProbaT,
																	discretizationBetweenDates) ;
	size_t nbSimu = 100000 ;
	mcQuadPricer->price(vanillaSwaption_Test_PTR, nbSimu) ;

}




VanillaSwaption_PTR setSwaptionATM_Quad(const CheyetteQuad_Model_PTR modele_test_PTR,
										const Tenor& floatTenor, const Tenor& fixedTenor, 
										const size_t a, const size_t b)
{
	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
	
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;

	double strike =  - 1000. ; //temporaire avant de mettre strike ATM
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	Cheyette_SwaptionPricer_QuadApprox approx(modele_test_PTR);

//calcul du strike ATM pour le swap
	double strikeATM = approx.swapRate0() ;
	swap.set_strike(strikeATM) ;
	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;

	return swaption_PTR_test_ATM ;
}